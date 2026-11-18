import test from 'node:test';
import assert from 'node:assert/strict';
import {
  G, TWO_PI, wavenumber, wavelength, deepWavelength, groupSpeed, groupFactor,
  shoalingCoefficient, stokes2Coefficient, windStressFactor, smbShallow, smbDeep,
  durationLimitedFetch, fetchGrowthTime, seichePeriod, equilibriumSetup, windStress,
  tmaFactor, frequencyBins, depthClass,
} from '../src/physics.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg}: ${a} vs ${b}`);

test('dispersion: deep-water limit k = omega^2/g', () => {
  const T = 8, w = TWO_PI / T;
  close(wavenumber(w, 5000), w * w / G, 1e-9, 'deep k');
  close(wavelength(T, 5000), deepWavelength(T), 1e-9, 'deep L');
});

test('dispersion: shallow-water limit k = omega / sqrt(g d)', () => {
  const T = 600, d = 10, w = TWO_PI / T;
  // exact kd·tanh(kd)=x² differs from the long-wave limit by O((kd)²/6)
  close(wavenumber(w, d), w / Math.sqrt(G * d), 1e-4, 'shallow k');
});

test('dispersion: solution satisfies omega^2 = g k tanh(kd) at intermediate depth', () => {
  for (const [T, d] of [[7.4, 19], [4, 3], [12, 40], [2, 1], [16, 250]]) {
    const w = TWO_PI / T;
    const k = wavenumber(w, d);
    close(G * k * Math.tanh(k * d), w * w, 1e-10, `residual T=${T} d=${d}`);
  }
});

test('group factor n goes from 1 (shallow) to 1/2 (deep)', () => {
  close(groupFactor(1e-4), 1, 1e-6, 'shallow n');
  close(groupFactor(30), 0.5, 1e-9, 'deep n');
});

test('shoaling: Ks dips to about 0.91 then grows past 1 toward shore', () => {
  const T = 10, w = TWO_PI / T;
  const L0 = deepWavelength(T);
  // Minimum of Ks is ~0.913 near d/L0 ≈ 0.16
  const ksMin = shoalingCoefficient(w, 5000, 0.157 * L0);
  close(ksMin, 0.913, 0.01, 'Ks minimum');
  const ksShallow = shoalingCoefficient(w, 5000, 0.02 * L0);
  close(ksShallow, 1.226, 0.01, 'Ks at d/L0 = 0.02 (table value 1.23)');
});

test('Stokes second-order coefficient → 1/2 in deep water', () => {
  close(stokes2Coefficient(10), 0.5, 1e-6, 'C2 deep');
  assert.ok(stokes2Coefficient(0.3) > 10, 'C2 grows in shallow water');
});

test('SMB shallow-water growth reproduces hand calculation for Lake Erie', () => {
  // U10 = 15 m/s, F = 300 km, d = 19 m  (worked by hand: Hs ≈ 3.10 m, T ≈ 7.38 s)
  const UA = windStressFactor(15);
  close(UA, 19.8, 0.01, 'UA');
  const { Hs, Tp } = smbShallow(UA, 300e3, 19);
  close(Hs, 3.10, 0.02, 'Hs');
  close(Tp, 7.38, 0.02, 'Tp');
  const deep = smbDeep(UA, 300e3);
  close(deep.Hs, 5.13, 0.02, 'deep Hs');
  close(deep.Tp, 9.59, 0.02, 'deep Tp');
});

test('SMB deep limit matches SPM deep-water coefficients for short fetch', () => {
  const UA = 20, F = 1000;
  const { Hs, Tp } = smbDeep(UA, F);
  const Fh = G * F / (UA * UA);
  close(Hs * G / (UA * UA), 1.6e-3 * Math.sqrt(Fh), 0.01, 'gH/UA^2');
  close(Tp * G / UA, 0.2857 * Math.cbrt(Fh), 0.01, 'gT/UA');
});

test('duration-limited fetch inverts fetchGrowthTime', () => {
  const UA = 19.8, F = 120e3;
  const t = fetchGrowthTime(UA, F);
  close(durationLimitedFetch(UA, t), F, 1e-9, 'round trip');
});

test('seiche period: uniform basin equals Merian formula', () => {
  const L = 388e3, d = 19, n = 400;
  const depths = new Float64Array(n).fill(d);
  const T1 = seichePeriod(depths, L / n);
  close(T1, 2 * L / Math.sqrt(G * d), 1e-9, 'Merian');
  close(T1 / 3600, 15.8, 0.01, 'about 15.8 h for an Erie-like channel');
});

test('equilibrium set-up: uniform basin tilt equals tau L / (rho g h)', () => {
  const L = 388e3, d = 19, n = 400;
  const depths = new Float64Array(n).fill(d);
  const tau = windStress(15);
  const eta = equilibriumSetup(depths, L / n, tau);
  const diff = eta[n - 1] - eta[0];
  close(diff, tau * L / (998 * G * d), 0.03, 'end-to-end set-up');
  close(diff, 1.0, 0.1, 'about 1 m for Lake Erie at 15 m/s');
  const mean = eta.reduce((a, b) => a + b, 0) / n;
  assert.ok(Math.abs(mean) < 1e-9, 'volume conserved');
});

test('TMA factor is 1 in deep water and small for long waves in shallow water', () => {
  assert.equal(tmaFactor(0.2, 500), 1);
  assert.ok(tmaFactor(0.05, 5) < 0.2);
});

test('frequency bins tile the range', () => {
  const { f, df } = frequencyBins(16, 0.05, 1.0);
  const total = df.reduce((a, b) => a + b, 0);
  const r = Math.pow(1.0 / 0.05, 1 / 15);
  close(total, 1.0 * Math.sqrt(r) - 0.05 / Math.sqrt(r), 1e-9, 'bin widths sum to range');
  assert.equal(f.length, 16);
});

test('depth classes', () => {
  assert.equal(depthClass(TWO_PI * 0.6), 'deep');
  assert.equal(depthClass(TWO_PI * 0.2), 'intermediate');
  assert.equal(depthClass(TWO_PI * 0.01), 'shallow');
});
