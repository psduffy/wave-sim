// slowmodel.js — one-dimensional shallow-water model of the whole basin.
//
// This is the "slow" physics: wind stress pushes the water toward the
// downwind shore (set-up); when the wind changes the tilted surface sloshes
// back as a seiche.  It is the depth-averaged, long-wave form of the
// momentum and continuity equations (Dean & Dalrymple ch. 5; Vreugdenhil,
// "Numerical Methods for Shallow-Water Flow" ch. 2):
//
//   ∂η/∂t = −∂(h u)/∂x
//   ∂u/∂t = −g ∂η/∂x + τ_wind/(ρ h) − C_f |u| u / h
//
// solved on a staggered (Arakawa C) grid with the forward–backward scheme:
// velocities first, then water levels with the updated velocities.  The
// seiche period falls out of the model itself; nothing is prescribed.

import { G, RHO_WATER } from './physics.js';

export const BOTTOM_DRAG = 0.0025;   // quadratic bottom-friction coefficient, typical for lakes
const H_MIN = 0.5;                   // m, keeps shore cells from drying out

export class SlowModel {
  /**
   * depths: still-water depth at N cell centres; dx: cell size (m);
   * widths: plan-view width at the same centres (m), or null for a uniform channel.
   * Width enters only the continuity equation: ∂η/∂t = −(1/b) ∂(b h u)/∂x.
   */
  constructor(depths, dx, widths = null) {
    this.n = depths.length;
    this.dx = dx;
    this.d = Float64Array.from(depths, (v) => Math.max(H_MIN, v));
    this.b = widths ? Float64Array.from(widths, (v) => Math.max(1, v)) : new Float64Array(this.n).fill(1);
    this.eta = new Float64Array(this.n);
    this.u = new Float64Array(this.n - 1);        // interior faces
    this.time = 0;                                // lake seconds
    let hmax = 0;
    for (const v of this.d) hmax = Math.max(hmax, v);
    this.dtMax = (0.6 * dx) / Math.sqrt(G * hmax);
    this.samples = [];                             // gauge history
    this.sampleEvery = 120;                        // s
    this.lastSample = -Infinity;
    this.historySeconds = 72 * 3600;
  }

  /** Replace the bed (same grid) without disturbing the water. */
  setDepths(depths, widths = null) {
    for (let j = 0; j < this.n; j++) this.d[j] = Math.max(H_MIN, depths[j]);
    if (widths) for (let j = 0; j < this.n; j++) this.b[j] = Math.max(1, widths[j]);
    let hmax = 0;
    for (const v of this.d) hmax = Math.max(hmax, v);
    this.dtMax = (0.6 * this.dx) / Math.sqrt(G * hmax);
  }

  /** Start from rest, or from a prescribed surface (e.g. equilibrium set-up). */
  reset(eta = null) {
    this.eta.fill(0);
    if (eta) for (let j = 0; j < this.n; j++) this.eta[j] = eta[j];
    this.u.fill(0);
    this.time = 0;
    this.samples.length = 0;
    this.lastSample = -Infinity;
    this.sample();
  }

  /** Advance by `dt` lake seconds under wind stress `tau` (Pa, + toward +x). */
  step(dt, tau) {
    const nSub = Math.max(1, Math.ceil(dt / this.dtMax));
    const h = dt / nSub;
    for (let s = 0; s < nSub; s++) this.substep(h, tau);
    this.time += dt;
    if (this.time - this.lastSample >= this.sampleEvery) this.sample();
  }

  substep(dt, tau) {
    const { n, dx, d, b, eta, u } = this;
    const g = G;
    const inv = 1 / (RHO_WATER);
    // momentum at faces
    for (let f = 0; f < n - 1; f++) {
      const hf = Math.max(H_MIN, 0.5 * (d[f] + d[f + 1] + eta[f] + eta[f + 1]));
      const uf = u[f];
      const accel = -g * (eta[f + 1] - eta[f]) / dx + (tau * inv) / hf - (BOTTOM_DRAG * Math.abs(uf) * uf) / hf;
      u[f] = uf + dt * accel;
    }
    // continuity at centres, closed walls at both ends
    let fluxLeft = 0;
    for (let j = 0; j < n; j++) {
      let fluxRight = 0;
      if (j < n - 1) {
        const hf = Math.max(H_MIN, 0.5 * (d[j] + d[j + 1] + eta[j] + eta[j + 1]));
        fluxRight = 0.5 * (b[j] + b[j + 1]) * hf * u[j];
      }
      eta[j] -= (dt / (dx * b[j])) * (fluxRight - fluxLeft);
      fluxLeft = fluxRight;
    }
  }

  sample() {
    this.lastSample = this.time;
    this.samples.push({ t: this.time, west: this.eta[0], east: this.eta[this.n - 1] });
    const cutoff = this.time - this.historySeconds;
    while (this.samples.length > 2 && this.samples[0].t < cutoff) this.samples.shift();
  }

  /** Water level at an arbitrary x (m), linear between cell centres. */
  levelAt(x) {
    const { n, dx, eta } = this;
    const p = x / dx - 0.5;
    if (p <= 0) return eta[0];
    if (p >= n - 1) return eta[n - 1];
    const j = Math.floor(p);
    const t = p - j;
    return eta[j] * (1 - t) + eta[j + 1] * t;
  }
}
