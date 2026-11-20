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
export const H_MIN = 0.5;            // m, the film left on a "dry" cell

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

  /**
   * Advance by `dt` lake seconds under wind stress `tau` (Pa, + toward +x).
   * The sub-step honours the CFL condition on the *current* state: long-wave
   * speed in the deepest water plus the fastest current, so a storm that
   * piles water 100 m deep against one shore stays stable.
   */
  step(dt, tau) {
    let remaining = dt;
    while (remaining > 1e-9) {
      let hmax = 0;
      for (let j = 0; j < this.n; j++) hmax = Math.max(hmax, this.d[j] + this.eta[j]);
      let umax = 0;
      for (let f = 0; f < this.n - 1; f++) umax = Math.max(umax, Math.abs(this.u[f]));
      const dtStable = (0.6 * this.dx) / (Math.sqrt(G * Math.max(H_MIN, hmax)) + umax);
      const h = Math.min(remaining, dtStable);
      this.substep(h, tau);
      remaining -= h;
    }
    this.time += dt;
    if (this.time - this.lastSample >= this.sampleEvery) this.sample();
  }

  substep(dt, tau) {
    const { n, dx, d, b, eta, u } = this;
    const g = G;
    const inv = 1 / RHO_WATER;
    const flux = this.flux || (this.flux = new Float64Array(n - 1));
    // momentum at faces; a face with no water on it carries no flow
    for (let f = 0; f < n - 1; f++) {
      const hf = 0.5 * (d[f] + d[f + 1] + eta[f] + eta[f + 1]);
      if (hf <= H_MIN * 1.001) { u[f] = 0; flux[f] = 0; continue; }
      const uf = u[f];
      const accel = -g * (eta[f + 1] - eta[f]) / dx + (tau * inv) / hf - (BOTTOM_DRAG * Math.abs(uf) * uf) / hf;
      u[f] = uf + dt * accel;
      flux[f] = 0.5 * (b[f] + b[f + 1]) * hf * u[f];
    }
    // a cell may not give away more water than it has above the dry limit:
    // scale its outgoing fluxes down together, which keeps the volume exact
    for (let j = 0; j < n; j++) {
      const out = (j < n - 1 && flux[j] > 0 ? flux[j] : 0) + (j > 0 && flux[j - 1] < 0 ? -flux[j - 1] : 0);
      if (out <= 0) continue;
      const available = Math.max(0, d[j] + eta[j] - H_MIN) * dx * b[j];
      if (out * dt > available) {
        const s = available / (out * dt);
        if (j < n - 1 && flux[j] > 0) { flux[j] *= s; u[j] *= s; }
        if (j > 0 && flux[j - 1] < 0) { flux[j - 1] *= s; u[j - 1] *= s; }
      }
    }
    // continuity at centres, closed walls at both ends
    for (let j = 0; j < n; j++) {
      const fluxRight = j < n - 1 ? flux[j] : 0;
      const fluxLeft = j > 0 ? flux[j - 1] : 0;
      eta[j] -= (dt / (dx * b[j])) * (fluxRight - fluxLeft);
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
