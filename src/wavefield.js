// wavefield.js — the wind-wave field along the section.
//
// Model in one sentence: the wind builds a sea whose height and period grow
// with fetch and are capped by the water depth it has blown over (SMB growth
// curves), and that sea is then carried over the local bed by linear wave
// theory (each spectral component keeps its frequency, changes its
// wavelength and speed with depth, shoals so that energy flux is conserved)
// until the depth-limited breaking criterion removes the excess.
//
// Simplifications, stated so they can be questioned:
//  • One horizontal dimension: no refraction, no directional spreading.
//  • The "generation depth" used in the growth curves is the mean depth
//    between the upwind shore and the point in question.
//  • Breaking is a hard cap on Hs (γ·d) with permanent energy loss, plus a
//    heuristic regrowth length so the wind can rebuild the sea past a sill.
//  • Components are summed linearly; the optional second-order crest
//    sharpening is applied per component and ignores cross terms.
//  • The phase of each component is integrated over the still-water depth.
//    Wind set-up changes the depth by centimetres, which would shift the
//    phase by hundreds of radians over a 400 km integration path if it were
//    included; a real wave train does not jump like that, its phase is
//    continuous in time and the local wavelength adjusts slowly.
//
// Performance: k, Cg and the Stokes coefficient of every component are
// tabulated against log(depth) once; the marches then only interpolate.

import {
  G, TWO_PI, wavenumber, groupFactor, jonswapShape, tmaFactor, frequencyBins,
  smbShallow, breakingLimitHs, stokes2Coefficient,
} from './physics.js';

const REGROW_LENGTH = 20e3;   // m, e-folding fetch for the sea to recover after breaking
const DRY_DEPTH = 0.3;        // m
const FILM_DEPTH = 1.0;       // m, water this thin is a blown-dry bed: no breaking marks
const STOKES_CAP = 0.3;       // max second-harmonic / first-harmonic amplitude ratio
const LUT_N = 256;
const LUT_MIN = Math.log(0.05);
const LUT_MAX = Math.log(3000);
const LUT_STEP = (LUT_MAX - LUT_MIN) / (LUT_N - 1);

export class WaveField {
  constructor({ nComp = 72, fMin = 0.02, fMax = 0.8, seed = 12345 } = {}) {
    this.nComp = nComp;
    const { f, df } = frequencyBins(nComp, fMin, fMax);
    this.f = f;
    this.df = df;
    this.omega = Float64Array.from(f, (v) => TWO_PI * v);
    this.phase0 = new Float64Array(nComp);
    let s = seed >>> 0;
    for (let i = 0; i < nComp; i++) {
      s = (1664525 * s + 1013904223) >>> 0;          // LCG, deterministic
      this.phase0[i] = (s / 4294967296) * TWO_PI;
    }
    this.buildLookup();
    this.coarse = null;
    this.fine = null;
  }

  /** k(d), Cg(d) and the Stokes coefficient per component on a log-depth grid. */
  buildLookup() {
    const N = this.nComp;
    this.kLut = new Float64Array(N * LUT_N);
    this.cgLut = new Float64Array(N * LUT_N);
    this.c2Lut = new Float64Array(N * LUT_N);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < LUT_N; j++) {
        const d = Math.exp(LUT_MIN + j * LUT_STEP);
        const k = wavenumber(this.omega[i], d);
        const idx = i * LUT_N + j;
        this.kLut[idx] = k;
        this.cgLut[idx] = (groupFactor(k * d) * this.omega[i]) / k;
        this.c2Lut[idx] = stokes2Coefficient(k * d);
      }
    }
  }

  /** Fractional table position for depth d (clamped to the table). */
  static depthIndex(d) {
    if (!(d > 0.05)) return 0;
    const p = (Math.log(d) - LUT_MIN) / LUT_STEP;
    return p >= LUT_N - 1 ? LUT_N - 1 - 1e-9 : p;
  }

  lookup(table, i, p) {
    const j = p | 0;
    const t = p - j;
    const base = i * LUT_N + j;
    return table[base] * (1 - t) + table[base + 1] * t;
  }

  kAt(i, p) { return this.lookup(this.kLut, i, p); }
  cgAt(i, p) { return this.lookup(this.cgLut, i, p); }
  c2At(i, p) { return this.lookup(this.c2Lut, i, p); }

  /**
   * Spectral weights (unnormalised energy per bin) for a sea of peak period
   * Tp at depth d: JONSWAP × TMA × bin width, written out with the
   * per-component powers precomputed (same result as jonswapShape/tmaFactor).
   */
  weights(Tp, dEff, out) {
    if (!this.fm4) {
      this.fm4 = Float64Array.from(this.f, (v) => Math.pow(v, -4));
      this.fm5df = Float64Array.from(this.f, (v, i) => Math.pow(v, -5) * this.df[i]);
    }
    const fp = 1 / Tp;
    const fp4 = fp * fp * fp * fp;
    const whScale = TWO_PI * Math.sqrt(dEff / G);
    const lnGamma = Math.log(3.3);
    let m0 = 0;
    for (let i = 0; i < this.nComp; i++) {
      const fi = this.f[i];
      const sigma = fi <= fp ? 0.07 : 0.09;
      const dfp = (fi - fp) / (sigma * fp);
      const r = Math.exp(-0.5 * dfp * dfp);
      let w = this.fm5df[i] * Math.exp(-1.25 * fp4 * this.fm4[i] + lnGamma * r);
      const wh = fi * whScale;
      if (wh <= 1) w *= 0.5 * wh * wh;
      else if (wh < 2) w *= 1 - 0.5 * (2 - wh) * (2 - wh);
      out[i] = w;
      m0 += w;
    }
    return m0;
  }

  /**
   * March the sea across the basin on the coarse grid.
   *   depths: still-water depth at M nodes, spacing dx
   *   setup:  mean water-level deviation at the same nodes (or null)
   *   UA:     effective wind-stress factor (m/s)
   *   dir:    +1 wind toward +x, −1 toward −x
   *   fetchCap: duration-limited fetch (m) — Infinity if fully developed
   *   opts:   { shoaling, breaking }
   * The phase table is rebuilt only when the bed, the grid or the direction
   * changed (it depends on nothing else).
   */
  computeCoarse(depths, dx, setup, UA, dir, fetchCap, opts) {
    const M = depths.length;
    const N = this.nComp;
    const c = this.coarse && this.coarse.M === M ? this.coarse : this.allocCoarse(M);
    const phaseDirty = c.dx !== dx || c.dir !== dir || !c.hasPhase || !sameDepths(c.still, depths);
    c.dx = dx;
    c.dir = dir;
    const { d, still, pStill, pTot, fetch, dEff, HsGen, Tp, Hs, phi, breaking, ampGen, Phi, Lp } = c;
    const w = c.scratch;

    for (let j = 0; j < M; j++) {
      still[j] = depths[j];
      d[j] = depths[j] + (setup ? setup[j] : 0);
      pStill[j] = WaveField.depthIndex(still[j]);
      pTot[j] = WaveField.depthIndex(d[j]);
    }

    let sumDepth = 0;
    let wet = 0;
    let phiMem = 1;
    const start = dir > 0 ? 0 : M - 1;
    const stop = dir > 0 ? M : -1;
    const upwindShore = dir > 0 ? 0 : (M - 1) * dx;
    for (let j = start; j !== stop; j += dir) {
      const dj = d[j];
      const isWet = dj > DRY_DEPTH;
      const F = Math.abs(j * dx - upwindShore);
      fetch[j] = F;
      if (isWet) { sumDepth += dj; wet++; }
      const de = wet > 0 ? sumDepth / wet : Math.max(DRY_DEPTH, dj);
      dEff[j] = de;
      const pEff = WaveField.depthIndex(de);

      const Feff = Math.min(F, fetchCap);
      const sea = smbShallow(UA, Feff, de);
      HsGen[j] = sea.Hs;
      Tp[j] = sea.Tp > 0 ? sea.Tp : 1;

      // generation spectrum → component amplitudes at the generation depth
      let m0 = 0;
      if (sea.Hs > 0 && isWet) m0 = this.weights(Tp[j], de, w);
      const scale = m0 > 0 ? (sea.Hs / 4) * Math.sqrt(2 / m0) : 0;
      for (let i = 0; i < N; i++) ampGen[i * M + j] = scale * Math.sqrt(w[i]);

      // transform to the local depth
      let e = 0;
      const pj = pTot[j];
      for (let i = 0; i < N; i++) {
        let a = ampGen[i * M + j];
        if (a > 0 && opts.shoaling) a *= Math.sqrt(this.cgAt(i, pEff) / this.cgAt(i, pj));
        a *= phiMem;
        w[i] = isWet ? a : 0;
        e += w[i] * w[i];
      }
      let hs = 4 * Math.sqrt(e / 2);
      const kp = isWet ? wavenumber(TWO_PI / Tp[j], dj) : Infinity;
      Lp[j] = isWet ? TWO_PI / kp : 0;
      let brk = 0;
      if (isWet && opts.breaking && hs > 0) {
        const lim = breakingLimitHs(dj, Lp[j], kp * dj);
        if (hs > lim) {
          const s = lim / hs;
          phiMem *= s;
          hs = lim;
          brk = dj > FILM_DEPTH ? 1 : 0;
        }
      }
      Hs[j] = hs;
      phi[j] = phiMem;
      breaking[j] = brk;
      phiMem += (1 - phiMem) * (1 - Math.exp(-dx / REGROW_LENGTH));
    }

    if (phaseDirty) {
      // cumulative ∫k dx in the propagation direction, over still water
      for (let j = start; j !== stop; j += dir) {
        const prev = j - dir;
        const wetHere = still[j] > DRY_DEPTH;
        if (j === start) { for (let i = 0; i < N; i++) Phi[i * M + j] = 0; continue; }
        const wetPrev = still[prev] > DRY_DEPTH;
        for (let i = 0; i < N; i++) {
          const kHere = wetHere ? this.kAt(i, pStill[j]) : 0;
          const kPrev = wetPrev ? this.kAt(i, pStill[prev]) : 0;
          Phi[i * M + j] = Phi[i * M + prev] + 0.5 * (kHere + kPrev) * dx;
        }
      }
      c.hasPhase = true;
      c.phaseVersion++;
    }
    this.coarse = c;
    return c;
  }

  allocCoarse(M) {
    const N = this.nComp;
    return {
      M, dx: 0, dir: 0, hasPhase: false, phaseVersion: 0,
      d: new Float64Array(M), still: new Float64Array(M), pStill: new Float64Array(M), pTot: new Float64Array(M),
      fetch: new Float64Array(M), dEff: new Float64Array(M),
      HsGen: new Float64Array(M), Tp: new Float64Array(M), Hs: new Float64Array(M),
      phi: new Float64Array(M), breaking: new Uint8Array(M), Lp: new Float64Array(M),
      ampGen: new Float64Array(N * M), Phi: new Float64Array(N * M),
      scratch: new Float64Array(N),
    };
  }

  /** Linear interpolation helper on the coarse grid. */
  coarseAt(arr, x) {
    const c = this.coarse;
    const p = x / c.dx;
    const j = Math.max(0, Math.min(c.M - 2, Math.floor(p)));
    const t = Math.max(0, Math.min(1, p - j));
    return arr[j] * (1 - t) + arr[j + 1] * t;
  }

  /**
   * Fine field for a window [x0, x1] sampled at n points.
   *   depthAt(x): still-water depth; setupAt(x): mean level (may be null)
   *   blend: keep the previous amplitude tables so surface() can fade from
   *          them to the new ones (false → the new tables replace both)
   * Returns per-point arrays and per-component amplitude/phase tables.
   */
  computeFine(x0, x1, n, depthAt, setupAt, opts, blend = true) {
    const c = this.coarse;
    const N = this.nComp;
    const f = this.fine && this.fine.n === n ? this.fine : this.allocFine(n);
    const sameWindow = f.x0 === x0 && f.x1 === x1 && f.phaseVersion === c.phaseVersion;
    const canBlend = blend && sameWindow && f.valid;
    // the tables just computed become "previous"; write the new ones into the other buffer
    [f.amp, f.ampPrev] = [f.ampPrev, f.amp];
    [f.amp2, f.amp2Prev] = [f.amp2Prev, f.amp2];
    [f.eta0, f.eta0Prev] = [f.eta0Prev, f.eta0];
    f.x0 = x0; f.x1 = x1;
    const { x, d, eta0, Hs, Tp, breaking, dEff, amp, amp2, Phi, k, pStill } = f;
    const dxf = (x1 - x0) / (n - 1);
    const dir = c.dir;
    const w = c.scratch;

    for (let p = 0; p < n; p++) {
      x[p] = x0 + p * dxf;
      eta0[p] = setupAt ? setupAt(x[p]) : 0;
      const still = depthAt(x[p]);
      d[p] = still + eta0[p];
      pStill[p] = WaveField.depthIndex(still);
      dEff[p] = this.coarseAt(c.dEff, x[p]);
      Tp[p] = this.coarseAt(c.Tp, x[p]);
    }

    // amplitudes: generation (interpolated) × shoaling × dissipation memory, then the breaking cap
    for (let p = 0; p < n; p++) {
      const dj = d[p];
      const isWet = dj > DRY_DEPTH;
      const phiMem = this.coarseAt(c.phi, x[p]);
      const pEff = WaveField.depthIndex(dEff[p]);
      const pj = WaveField.depthIndex(dj);
      const pc = x[p] / c.dx;
      const j0 = Math.max(0, Math.min(c.M - 2, Math.floor(pc)));
      const t = Math.max(0, Math.min(1, pc - j0));
      let e = 0;
      for (let i = 0; i < N; i++) {
        let a = c.ampGen[i * c.M + j0] * (1 - t) + c.ampGen[i * c.M + j0 + 1] * t;
        if (isWet) {
          if (a > 0 && opts.shoaling) a *= Math.sqrt(this.cgAt(i, pEff) / this.cgAt(i, pj));
          a *= phiMem;
        } else a = 0;
        w[i] = a;
        e += a * a;
      }
      let hs = 4 * Math.sqrt(e / 2);
      let brk = 0;
      if (isWet && opts.breaking && hs > 0) {
        const kp = wavenumber(TWO_PI / Tp[p], dj);
        const lim = breakingLimitHs(dj, TWO_PI / kp, kp * dj);
        if (hs > lim) {
          const s = lim / hs;
          for (let i = 0; i < N; i++) w[i] *= s;
          hs = lim;
          brk = dj > FILM_DEPTH ? 1 : 0;
        }
      }
      Hs[p] = hs;
      breaking[p] = brk;
      for (let i = 0; i < N; i++) {
        const idx = i * n + p;
        const ki = isWet ? this.kAt(i, pj) : 0;
        k[idx] = ki;
        // a component shorter than ~2.5 samples cannot be drawn, only aliased: leave it out of the picture
        const drawable = ki * dxf <= 2.5;
        amp[idx] = drawable ? w[i] : 0;
        let a2 = 0;
        if (opts.nonlinear && drawable && w[i] > 0 && isWet) {
          const ratio = Math.min(STOKES_CAP, ki * w[i] * this.c2At(i, pj));
          a2 = ratio * w[i];
        }
        amp2[idx] = a2;
      }
    }

    // phase: from the nearest coarse node on the upwind side, over still water
    if (!sameWindow) {
      const jRef = dir > 0 ? Math.max(0, Math.floor(x0 / c.dx)) : Math.min(c.M - 1, Math.ceil(x1 / c.dx));
      const xRef = jRef * c.dx;
      const pStart = dir > 0 ? 0 : n - 1;
      const gap = Math.abs(x[pStart] - xRef);
      const nGap = 24;
      const gapP = new Float64Array(nGap + 1);
      for (let s = 0; s <= nGap; s++) gapP[s] = WaveField.depthIndex(depthAt(xRef + dir * (gap * s) / nGap));
      for (let i = 0; i < N; i++) {
        let acc = c.Phi[i * c.M + jRef];
        let kPrev = this.kAt(i, gapP[0]);
        for (let s = 1; s <= nGap; s++) {
          const ks = this.kAt(i, gapP[s]);
          acc += 0.5 * (kPrev + ks) * (gap / nGap);
          kPrev = ks;
        }
        Phi[i * n + pStart] = acc;
        for (let p = pStart + dir; p >= 0 && p < n; p += dir) {
          const prev = p - dir;
          Phi[i * n + p] = Phi[i * n + prev] + 0.5 * (this.kAt(i, pStill[p]) + this.kAt(i, pStill[prev])) * dxf;
        }
      }
      f.phaseVersion = c.phaseVersion;
    }

    if (!canBlend) {
      f.ampPrev.set(amp);
      f.amp2Prev.set(amp2);
      f.eta0Prev.set(eta0);
    }
    f.valid = true;
    this.fine = f;
    return f;
  }

  allocFine(n) {
    const N = this.nComp;
    return {
      n, x0: NaN, x1: NaN, valid: false, phaseVersion: -1,
      x: new Float64Array(n), d: new Float64Array(n), pStill: new Float64Array(n),
      eta0: new Float64Array(n), eta0Prev: new Float64Array(n),
      Hs: new Float64Array(n), Tp: new Float64Array(n), breaking: new Uint8Array(n),
      dEff: new Float64Array(n),
      amp: new Float64Array(N * n), ampPrev: new Float64Array(N * n),
      amp2: new Float64Array(N * n), amp2Prev: new Float64Array(N * n),
      Phi: new Float64Array(N * n), k: new Float64Array(N * n),
      eta: new Float64Array(n),
    };
  }

  /**
   * Surface elevation (relative to still water, including set-up) at time t.
   * blendW in [0, 1] fades from the previous amplitude tables (0) to the
   * current ones (1), so a rebuild never shows as a step.
   */
  surface(t, blendW = 1) {
    const f = this.fine;
    const { n, amp, ampPrev, amp2, amp2Prev, Phi, eta, eta0, eta0Prev } = f;
    const N = this.nComp;
    const wB = Math.max(0, Math.min(1, blendW));
    const wA = 1 - wB;
    for (let p = 0; p < n; p++) eta[p] = eta0Prev[p] * wA + eta0[p] * wB;
    for (let i = 0; i < N; i++) {
      const wt = this.omega[i] * t - this.phase0[i];
      const base = i * n;
      for (let p = 0; p < n; p++) {
        const a = ampPrev[base + p] * wA + amp[base + p] * wB;
        if (a === 0) continue;
        const cs = Math.cos(Phi[base + p] - wt);
        eta[p] += a * cs + (amp2Prev[base + p] * wA + amp2[base + p] * wB) * (2 * cs * cs - 1);
      }
    }
    return eta;
  }

  /**
   * Orbital displacement of a water particle whose rest position is (x, z)
   * (z negative down from still water) at time t — sum over components.
   */
  particleDisplacement(p, z, t, blendW = 1) {
    const f = this.fine;
    const { n, amp, ampPrev, Phi, k, d, eta0, eta0Prev } = f;
    const N = this.nComp;
    const wB = Math.max(0, Math.min(1, blendW));
    const wA = 1 - wB;
    let dx = 0;
    let dz = 0;
    const depth = d[p];
    // z is measured from the still-water datum; the orbital profile wants the
    // elevation below the local mean surface
    const zz = Math.max(-depth, Math.min(0, z - (eta0Prev[p] * wA + eta0[p] * wB)));
    for (let i = 0; i < N; i++) {
      const idx = i * n + p;
      const a = ampPrev[idx] * wA + amp[idx] * wB;
      if (a === 0) continue;
      const ki = k[idx];
      const kd = ki * depth;
      let A;
      let B;
      if (kd > 20) {
        const e = Math.exp(ki * zz);
        A = a * e; B = a * e;
      } else {
        const s = Math.sinh(kd);
        A = (a * Math.cosh(ki * (zz + depth))) / s;
        B = (a * Math.sinh(ki * (zz + depth))) / s;
      }
      const th = Phi[idx] - this.omega[i] * t + this.phase0[i];
      dx += -A * Math.sin(th) * this.coarse.dir;
      dz += B * Math.cos(th);
    }
    return { dx, dz };
  }
}

function sameDepths(a, b) {
  if (a.length !== b.length) return false;
  for (let j = 0; j < a.length; j++) if (a[j] !== b[j]) return false;
  return true;
}
