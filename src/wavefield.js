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

import {
  G, TWO_PI, wavenumber, groupFactor, jonswapShape, tmaFactor, frequencyBins,
  smbShallow, breakingLimitHs, stokes2Coefficient,
} from './physics.js';

const REGROW_LENGTH = 20e3;   // m, e-folding fetch for the sea to recover after breaking
const DRY_DEPTH = 0.3;        // m
const STOKES_CAP = 0.3;       // max second-harmonic / first-harmonic amplitude ratio

export class WaveField {
  constructor({ nComp = 64, fMin = 0.04, fMax = 0.8, seed = 12345 } = {}) {
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
  }

  /** k(d) lookup per component on a log-depth grid; interpolated in log d. */
  buildLookup() {
    const nd = 240;
    this.lutMin = Math.log(0.05);
    this.lutMax = Math.log(3000);
    this.lutStep = (this.lutMax - this.lutMin) / (nd - 1);
    this.lut = new Float64Array(this.nComp * nd);
    this.lutN = nd;
    for (let i = 0; i < this.nComp; i++) {
      for (let j = 0; j < nd; j++) {
        const d = Math.exp(this.lutMin + j * this.lutStep);
        this.lut[i * nd + j] = wavenumber(this.omega[i], d);
      }
    }
  }

  kOf(i, d) {
    if (d <= 0.05) d = 0.05;
    if (d >= 3000) return (this.omega[i] * this.omega[i]) / G;
    const p = (Math.log(d) - this.lutMin) / this.lutStep;
    const j = Math.floor(p);
    const t = p - j;
    const base = i * this.lutN;
    return this.lut[base + j] * (1 - t) + this.lut[base + j + 1] * t;
  }

  cgOf(i, d) {
    const k = this.kOf(i, d);
    return (groupFactor(k * d) * this.omega[i]) / k;
  }

  /** Spectral weights (unnormalised energy per bin) for a sea of peak period Tp at depth d. */
  weights(Tp, dEff, out) {
    const fp = 1 / Tp;
    let m0 = 0;
    for (let i = 0; i < this.nComp; i++) {
      const w = jonswapShape(this.f[i], fp) * tmaFactor(this.f[i], dEff) * this.df[i];
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
   */
  computeCoarse(depths, dx, setup, UA, dir, fetchCap, opts) {
    const M = depths.length;
    const N = this.nComp;
    const c = this.coarse && this.coarse.M === M ? this.coarse : this.allocCoarse(M);
    c.dx = dx;
    c.dir = dir;
    const { d, fetch, dEff, HsGen, Tp, Hs, phi, breaking, ampGen, Phi, Lp } = c;
    const w = c.scratch;

    for (let j = 0; j < M; j++) d[j] = depths[j] + (setup ? setup[j] : 0);

    let sumDepth = 0;
    let wet = 0;
    let phiMem = 1;
    const start = dir > 0 ? 0 : M - 1;
    const stop = dir > 0 ? M : -1;
    for (let j = start; j !== stop; j += dir) {
      const dj = d[j];
      const isWet = dj > DRY_DEPTH;
      const upwindShore = dir > 0 ? 0 : (M - 1) * dx;
      const F = Math.abs(j * dx - upwindShore);
      fetch[j] = F;
      if (isWet) { sumDepth += dj; wet++; }
      const de = wet > 0 ? sumDepth / wet : Math.max(DRY_DEPTH, dj);
      dEff[j] = de;

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
      for (let i = 0; i < N; i++) {
        let a = ampGen[i * M + j];
        if (a > 0 && opts.shoaling) a *= Math.sqrt(this.cgOf(i, de) / this.cgOf(i, dj));
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
          brk = 1;
        }
      }
      Hs[j] = hs;
      phi[j] = phiMem;
      breaking[j] = brk;
      phiMem += (1 - phiMem) * (1 - Math.exp(-dx / REGROW_LENGTH));

      // phase: cumulative ∫k dx in the propagation direction
      const prev = j - dir;
      for (let i = 0; i < N; i++) {
        const kHere = isWet ? this.kOf(i, dj) : 0;
        if (prev === start - dir) {
          Phi[i * M + j] = 0;
        } else {
          const dprev = d[prev];
          const kPrev = dprev > DRY_DEPTH ? this.kOf(i, dprev) : 0;
          Phi[i * M + j] = Phi[i * M + prev] + 0.5 * (kHere + kPrev) * dx;
        }
      }
    }
    this.coarse = c;
    return c;
  }

  allocCoarse(M) {
    const N = this.nComp;
    return {
      M,
      d: new Float64Array(M), fetch: new Float64Array(M), dEff: new Float64Array(M),
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
   * Returns per-point arrays and per-component amplitude/phase tables.
   */
  computeFine(x0, x1, n, depthAt, setupAt, opts) {
    const c = this.coarse;
    const N = this.nComp;
    const f = this.fine && this.fine.n === n ? this.fine : this.allocFine(n);
    f.x0 = x0; f.x1 = x1;
    const { x, d, eta0, Hs, Tp, breaking, dEff, amp, amp2, Phi, k } = f;
    const dxf = (x1 - x0) / (n - 1);
    const dir = c.dir;
    const w = c.scratch;

    for (let p = 0; p < n; p++) {
      x[p] = x0 + p * dxf;
      eta0[p] = setupAt ? setupAt(x[p]) : 0;
      d[p] = depthAt(x[p]) + eta0[p];
      dEff[p] = this.coarseAt(c.dEff, x[p]);
      Tp[p] = this.coarseAt(c.Tp, x[p]);
    }

    // amplitudes: generation (interpolated) × shoaling × dissipation memory, then the breaking cap
    for (let p = 0; p < n; p++) {
      const dj = d[p];
      const isWet = dj > DRY_DEPTH;
      const phiMem = this.coarseAt(c.phi, x[p]);
      const de = dEff[p];
      const pc = x[p] / c.dx;
      const j0 = Math.max(0, Math.min(c.M - 2, Math.floor(pc)));
      const t = Math.max(0, Math.min(1, pc - j0));
      let e = 0;
      for (let i = 0; i < N; i++) {
        let a = c.ampGen[i * c.M + j0] * (1 - t) + c.ampGen[i * c.M + j0 + 1] * t;
        if (isWet) {
          if (a > 0 && opts.shoaling) a *= Math.sqrt(this.cgOf(i, de) / this.cgOf(i, dj));
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
          brk = 1;
        }
      }
      Hs[p] = hs;
      breaking[p] = brk;
      for (let i = 0; i < N; i++) {
        const idx = i * n + p;
        const ki = isWet ? this.kOf(i, dj) : 0;
        k[idx] = ki;
        amp[idx] = w[i];
        let a2 = 0;
        if (opts.nonlinear && w[i] > 0 && isWet) {
          const ratio = Math.min(STOKES_CAP, ki * w[i] * stokes2Coefficient(ki * dj));
          a2 = ratio * w[i];
        }
        amp2[idx] = a2;
      }
    }

    // phase: start from the nearest coarse node on the upwind side of the window
    const jRef = dir > 0 ? Math.max(0, Math.floor(x0 / c.dx)) : Math.min(c.M - 1, Math.ceil(x1 / c.dx));
    const xRef = jRef * c.dx;
    const pStart = dir > 0 ? 0 : n - 1;
    const gap = Math.abs(x[pStart] - xRef);
    const nGap = 24;
    for (let i = 0; i < N; i++) {
      // bridge the gap between the reference node and the window edge
      let acc = c.Phi[i * c.M + jRef];
      let kPrev = this.kOf(i, Math.max(DRY_DEPTH, c.d[jRef]));
      for (let s = 1; s <= nGap; s++) {
        const xs = xRef + dir * (gap * s) / nGap;
        const ds = depthAt(xs) + (setupAt ? setupAt(xs) : 0);
        const ks = ds > DRY_DEPTH ? this.kOf(i, ds) : 0;
        acc += 0.5 * (kPrev + ks) * (gap / nGap);
        kPrev = ks;
      }
      Phi[i * n + pStart] = acc;
      for (let p = pStart + dir; p >= 0 && p < n; p += dir) {
        const prev = p - dir;
        Phi[i * n + p] = Phi[i * n + prev] + 0.5 * (k[i * n + p] + k[i * n + prev]) * dxf;
      }
    }
    this.fine = f;
    return f;
  }

  allocFine(n) {
    const N = this.nComp;
    return {
      n,
      x: new Float64Array(n), d: new Float64Array(n), eta0: new Float64Array(n),
      Hs: new Float64Array(n), Tp: new Float64Array(n), breaking: new Uint8Array(n),
      dEff: new Float64Array(n),
      amp: new Float64Array(N * n), amp2: new Float64Array(N * n),
      Phi: new Float64Array(N * n), k: new Float64Array(N * n),
      eta: new Float64Array(n),
    };
  }

  /** Surface elevation (relative to still water, including set-up) at time t. */
  surface(t) {
    const f = this.fine;
    const { n, amp, amp2, Phi, eta, eta0 } = f;
    const N = this.nComp;
    for (let p = 0; p < n; p++) eta[p] = eta0[p];
    for (let i = 0; i < N; i++) {
      const wt = this.omega[i] * t - this.phase0[i];
      const base = i * n;
      for (let p = 0; p < n; p++) {
        const a = amp[base + p];
        if (a === 0) continue;
        const cs = Math.cos(Phi[base + p] - wt);
        eta[p] += a * cs + amp2[base + p] * (2 * cs * cs - 1);
      }
    }
    return eta;
  }

  /**
   * Orbital displacement of a water particle whose rest position is (x, z)
   * (z negative down from still water) at time t — sum over components.
   */
  particleDisplacement(p, z, t) {
    const f = this.fine;
    const { n, amp, Phi, k, d } = f;
    const N = this.nComp;
    let dx = 0;
    let dz = 0;
    const depth = d[p];
    const zz = Math.max(-depth, Math.min(0, z));
    for (let i = 0; i < N; i++) {
      const idx = i * n + p;
      const a = amp[idx];
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
