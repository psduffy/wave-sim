// physics.js — pure wave physics. No DOM, no state.
//
// Every formula here is a textbook or manual result; the source is cited next
// to it so the chain of reasoning can be followed down one more layer.
//
//   [AIRY]   Linear (Airy) wave theory: Dean & Dalrymple, "Water Wave Mechanics
//            for Engineers and Scientists" (1991), ch. 3–4.
//   [SPM]    U.S. Army Corps of Engineers, "Shore Protection Manual" (1984),
//            vol. I, ch. 3 (wind-wave growth, eqs. 3-28a, 3-33..3-35, 3-39..3-41).
//   [JONSWAP] Hasselmann et al. (1973); TMA depth factor: Bouws et al. (1985).
//   [WU]     Wu (1982), J. Geophys. Res. 87: wind-stress drag coefficient.
//   [MICHE]  Miche (1944) steepness limit; depth-limited breaking after
//            McCowan (1894) and the Coastal Engineering Manual (2002) part II-4.
//   [MERIAN] Merian (1828) / Lamb, "Hydrodynamics" §§ 183–185: seiche period.

export const G = 9.81;          // m/s^2
export const RHO_WATER = 998;   // kg/m^3, fresh water near 15 °C
export const RHO_AIR = 1.225;   // kg/m^3, sea-level air
export const TWO_PI = 2 * Math.PI;

// ---------------------------------------------------------------------------
// Dispersion: omega^2 = g k tanh(k d)   [AIRY]
// ---------------------------------------------------------------------------

/**
 * Wavenumber k (rad/m) for angular frequency omega (rad/s) at depth d (m).
 * Starts from Guo's (2002) explicit approximation and polishes with Newton
 * steps on f(kd) = kd·tanh(kd) − omega²d/g.
 */
export function wavenumber(omega, d) {
  if (!(d > 0)) return Infinity;
  if (omega <= 0) return 0;
  const x = omega * Math.sqrt(d / G);       // sqrt(omega^2 d / g)
  const target = x * x;
  let kd = target / Math.pow(1 - Math.exp(-Math.pow(x, 2.5)), 0.4);
  if (!isFinite(kd) || kd <= 0) kd = Math.max(x, target);
  for (let i = 0; i < 4; i++) {
    const th = Math.tanh(kd);
    const f = kd * th - target;
    const df = th + kd * (1 - th * th);
    const step = f / df;
    kd -= step;
    if (Math.abs(step) < 1e-12 * kd) break;
  }
  return kd / d;
}

/** n = Cg / C = ½ (1 + 2kd / sinh 2kd)   [AIRY] */
export function groupFactor(kd) {
  if (kd > 20) return 0.5;
  if (kd < 1e-6) return 1;
  return 0.5 * (1 + (2 * kd) / Math.sinh(2 * kd));
}

export function phaseSpeed(omega, k) { return omega / k; }
export function groupSpeed(omega, k, d) { return groupFactor(k * d) * omega / k; }

export function deepWavelength(T) { return (G * T * T) / TWO_PI; }
export function wavelength(T, d) { return TWO_PI / wavenumber(TWO_PI / T, d); }

/** Depth class by kd (d/L): deep if d/L > ½, shallow if d/L < 1/20.  [AIRY] */
export function depthClass(kd) {
  const dOverL = kd / TWO_PI;
  if (dOverL > 0.5) return 'deep';
  if (dOverL < 0.05) return 'shallow';
  return 'intermediate';
}

/**
 * Shoaling coefficient Ks = sqrt(Cg_from / Cg_to): energy flux E·Cg is
 * conserved along the ray, and E ∝ H², so H grows where Cg shrinks.  [AIRY]
 */
export function shoalingCoefficient(omega, dFrom, dTo) {
  const k1 = wavenumber(omega, dFrom);
  const k2 = wavenumber(omega, dTo);
  return Math.sqrt(groupSpeed(omega, k1, dFrom) / groupSpeed(omega, k2, dTo));
}

/**
 * Second-order Stokes crest coefficient: eta2 = k a² · C2(kd) · cos 2θ.
 * C2 → ½ in deep water and grows without bound as kd → 0 (where Stokes
 * theory stops being valid; see ursell()).  Dean & Dalrymple eq. 11.37.
 */
export function stokes2Coefficient(kd) {
  if (kd > 15) return 0.5;
  const s = Math.sinh(kd);
  const c = Math.cosh(kd);
  return (c * (2 + Math.cosh(2 * kd))) / (4 * s * s * s);
}

/** Ursell number H L² / d³; Stokes theory is adequate below ~26. */
export function ursell(H, L, d) { return (H * L * L) / (d * d * d); }

/**
 * Linear-theory orbital semi-axes at elevation z (negative down, 0 = still
 * water level) for a wave of amplitude a: horizontal A, vertical B.  [AIRY]
 */
export function orbitalSemiAxes(a, k, d, z) {
  const kd = k * d;
  const zz = Math.max(-d, Math.min(0, z));
  if (kd > 20) {
    const e = Math.exp(k * zz);
    return { A: a * e, B: a * e };
  }
  const s = Math.sinh(kd);
  return {
    A: (a * Math.cosh(k * (zz + d))) / s,
    B: (a * Math.sinh(k * (zz + d))) / s,
  };
}

// ---------------------------------------------------------------------------
// Wind
// ---------------------------------------------------------------------------

/** Wind-stress factor U_A = 0.71 U10^1.23 (U10 in m/s).  [SPM eq. 3-28a] */
export function windStressFactor(U10) { return U10 > 0 ? 0.71 * Math.pow(U10, 1.23) : 0; }

/** Drag coefficient C_D = (0.8 + 0.065 U10)·10⁻³.  [WU] */
export function dragCoefficient(U10) { return (0.8 + 0.065 * U10) * 1e-3; }

/** Surface wind stress τ = ρ_air C_D U10² (Pa). */
export function windStress(U10) { return RHO_AIR * dragCoefficient(U10) * U10 * U10; }

export function beaufort(U10) {
  const limits = [0.3, 1.6, 3.4, 5.5, 8.0, 10.8, 13.9, 17.2, 20.8, 24.5, 28.5, 32.7];
  let b = 0;
  while (b < limits.length && U10 >= limits[b]) b++;
  return b;
}

// ---------------------------------------------------------------------------
// Fetch- and depth-limited wind-wave growth (SMB method)  [SPM eqs. 3-39, 3-40]
// ---------------------------------------------------------------------------

/**
 * Significant wave height Hs (m) and period T (s) for wind-stress factor UA
 * (m/s), fetch F (m) and water depth d (m).  Deep water is the limit d → ∞
 * (both tanh depth factors → 1), which reproduces SPM eqs. 3-33/3-34.
 */
export function smbShallow(UA, F, d) {
  if (!(UA > 0) || !(F > 0) || !(d > 0)) return { Hs: 0, Tp: 0 };
  const UA2 = UA * UA;
  const Fh = (G * F) / UA2;         // dimensionless fetch
  const dh = (G * d) / UA2;         // dimensionless depth
  const A = Math.tanh(0.530 * Math.pow(dh, 0.75));
  const B = Math.tanh(0.833 * Math.pow(dh, 0.375));
  const Hh = 0.283 * A * Math.tanh((0.00565 * Math.sqrt(Fh)) / A);
  const Th = 7.54 * B * Math.tanh((0.0379 * Math.cbrt(Fh)) / B);
  return { Hs: (Hh * UA2) / G, Tp: (Th * UA) / G };
}

/** Deep-water counterpart (d → ∞) for comparisons. */
export function smbDeep(UA, F) { return smbShallow(UA, F, 1e6); }

/**
 * Duration limit: the fetch that a wind of stress factor UA can "fill" in
 * time t (s).  gt/UA = 68.8 (gF/UA²)^(2/3), inverted for F.  [SPM eq. 3-35]
 */
export function durationLimitedFetch(UA, t) {
  if (!(UA > 0) || !(t > 0)) return 0;
  return (Math.pow((G * t) / (68.8 * UA), 1.5) * UA * UA) / G;
}

/** Time (s) needed for the sea to become fetch-limited at fetch F.  [SPM eq. 3-35] */
export function fetchGrowthTime(UA, F) {
  if (!(UA > 0) || !(F > 0)) return 0;
  return 68.8 * (UA / G) * Math.pow((G * F) / (UA * UA), 2 / 3);
}

// ---------------------------------------------------------------------------
// Spectrum shape  [JONSWAP] with TMA depth factor
// ---------------------------------------------------------------------------

/** Unnormalised JONSWAP spectral shape at frequency f (Hz) for peak fp. */
export function jonswapShape(f, fp, gamma = 3.3) {
  if (!(f > 0) || !(fp > 0)) return 0;
  const sigma = f <= fp ? 0.07 : 0.09;
  const r = Math.exp(-((f - fp) * (f - fp)) / (2 * sigma * sigma * fp * fp));
  return Math.pow(f, -5) * Math.exp(-1.25 * Math.pow(f / fp, -4)) * Math.pow(gamma, r);
}

/**
 * Kitaigorodskii depth function φ(ω_h), ω_h = 2πf·sqrt(d/g).  Finite depth
 * removes energy from the low-frequency (long-wave) end of the spectrum.
 */
export function tmaFactor(f, d) {
  const wh = TWO_PI * f * Math.sqrt(d / G);
  if (wh <= 1) return 0.5 * wh * wh;
  if (wh < 2) return 1 - 0.5 * (2 - wh) * (2 - wh);
  return 1;
}

/** Log-spaced frequency bins from fMin to fMax (Hz): { f, df } arrays. */
export function frequencyBins(n, fMin, fMax) {
  const f = new Float64Array(n);
  const df = new Float64Array(n);
  const r = Math.pow(fMax / fMin, 1 / (n - 1));
  for (let i = 0; i < n; i++) f[i] = fMin * Math.pow(r, i);
  for (let i = 0; i < n; i++) {
    const lo = i === 0 ? f[0] / Math.sqrt(r) : Math.sqrt(f[i - 1] * f[i]);
    const hi = i === n - 1 ? f[n - 1] * Math.sqrt(r) : Math.sqrt(f[i] * f[i + 1]);
    df[i] = hi - lo;
  }
  return { f, df };
}

// ---------------------------------------------------------------------------
// Breaking
// ---------------------------------------------------------------------------

/** Depth-limited breaking index for Hs (individual waves break near H/d ≈ 0.78). */
export const GAMMA_SIG = 0.6;
/** Miche steepness limit applied to Hs with the peak wavelength. */
export const STEEPNESS_LIMIT = 0.10;

export function breakingLimitHs(d, Lp, kpd) {
  const byDepth = GAMMA_SIG * d;
  const bySteepness = STEEPNESS_LIMIT * Lp * Math.tanh(kpd);
  return Math.min(byDepth, bySteepness);
}

// ---------------------------------------------------------------------------
// Basin-scale response
// ---------------------------------------------------------------------------

/**
 * Equilibrium wind set-up: dη/dx = τ / (ρ g h), integrated along the basin
 * and shifted to conserve volume.  depths: array of still-water depths (m)
 * at spacing dx; tau: stress (Pa) positive toward +x.  Dean & Dalrymple §5.
 */
export function equilibriumSetup(depths, dx, tau, iterations = 3, widths = null) {
  const n = depths.length;
  const eta = new Float64Array(n);
  for (let it = 0; it < iterations; it++) {
    let acc = 0;
    const next = new Float64Array(n);
    for (let j = 1; j < n; j++) {
      const h = Math.max(0.5, 0.5 * (depths[j] + depths[j - 1] + eta[j] + eta[j - 1]));
      acc += (tau / (RHO_WATER * G * h)) * dx;
      next[j] = acc;
    }
    // shift so that the displaced volume Σ b·η·dx is zero
    let num = 0;
    let den = 0;
    for (let j = 0; j < n; j++) {
      const b = widths ? widths[j] : 1;
      num += b * next[j];
      den += b;
    }
    const mean = num / den;
    for (let j = 0; j < n; j++) eta[j] = next[j] - mean;
  }
  return eta;
}

/**
 * Fundamental seiche period T1 = 2 ∫ dx / sqrt(g h): the time for a long wave
 * to cross the basin and back.  Equals Merian's 2L/sqrt(gh) for uniform depth.
 */
export function seichePeriod(depths, dx) {
  let t = 0;
  for (let j = 0; j < depths.length; j++) {
    const h = Math.max(0.3, depths[j]);
    t += dx / Math.sqrt(G * h);
  }
  return 2 * t;
}

export function meanDepth(depths) {
  let s = 0;
  for (let j = 0; j < depths.length; j++) s += depths[j];
  return s / depths.length;
}
