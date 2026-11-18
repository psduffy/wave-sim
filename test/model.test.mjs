import test from 'node:test';
import assert from 'node:assert/strict';
import { PRESETS, presetById, sampleDepths, depthFunction, pchip } from '../src/bathymetry.js';
import { WaveField } from '../src/wavefield.js';
import { SlowModel } from '../src/slowmodel.js';
import { windStressFactor, windStress, equilibriumSetup, seichePeriod, G } from '../src/physics.js';

test('pchip passes through control points and never overshoots', () => {
  const f = pchip([0, 1, 2, 3], [0, 10, 10, 0]);
  assert.equal(f(1), 10);
  assert.equal(f(2), 10);
  for (let x = 0; x <= 3; x += 0.01) assert.ok(f(x) <= 10 + 1e-9 && f(x) >= -1e-9, `overshoot at ${x}`);
});

test('every preset starts and ends at the shore and has a sensible max depth', () => {
  for (const p of PRESETS) {
    const d = sampleDepths(p, 400);
    assert.equal(d[0], 0, p.id);
    assert.equal(d[d.length - 1], 0, p.id);
    const max = Math.max(...d);
    assert.ok(max > 5 && max < 500, `${p.id} max ${max}`);
  }
});

test('wave field: Erie, 15 m/s westerly — grows along the fetch, breaks only near the east shore', () => {
  const erie = presetById('erie');
  const M = 600;
  const depths = sampleDepths(erie, M);
  const dx = erie.length / (M - 1);
  const wf = new WaveField();
  const UA = windStressFactor(15);
  const c = wf.computeCoarse(depths, dx, null, UA, +1, Infinity, { shoaling: true, breaking: true });
  // monotone-ish growth over the first 200 km
  const at = (km) => c.Hs[Math.round((km * 1e3) / dx)];
  assert.ok(at(20) < at(100) && at(100) < at(200), `growth ${at(20)} ${at(100)} ${at(200)}`);
  // open-water height in the right range for a strong Erie westerly
  assert.ok(at(300) > 2.5 && at(300) < 4.0, `Hs at 300 km = ${at(300)}`);
  // breaking confined to the last few km
  let firstBreak = -1;
  for (let j = 0; j < M; j++) if (c.breaking[j]) { firstBreak = j; break; }
  assert.ok(firstBreak > 0, 'there is a breaking zone');
  assert.ok((M - 1 - firstBreak) * dx < 8e3, `breaking starts ${(M - 1 - firstBreak) * dx / 1e3} km from shore`);
  // no breaking in the open lake
  for (let j = 50; j < M - 20; j++) assert.equal(c.breaking[j], 0, `unexpected breaking at ${j * dx / 1e3} km`);
});

test('wave field: deep flat basin gives higher, longer waves than the shallow one', () => {
  const run = (id) => {
    const p = presetById(id);
    const M = 400;
    const depths = sampleDepths(p, M);
    const dx = p.length / (M - 1);
    const wf = new WaveField();
    const c = wf.computeCoarse(depths, dx, null, windStressFactor(15), +1, Infinity, { shoaling: true, breaking: true });
    const j = Math.round(180e3 / dx);
    return { Hs: c.Hs[j], Tp: c.Tp[j], Lp: c.Lp[j] };
  };
  const deep = run('flat-deep');
  const shallow = run('flat-shallow');
  assert.ok(deep.Hs > shallow.Hs, `Hs deep ${deep.Hs} shallow ${shallow.Hs}`);
  assert.ok(deep.Tp > shallow.Tp, `Tp deep ${deep.Tp} shallow ${shallow.Tp}`);
  assert.ok(deep.Hs / deep.Lp < shallow.Hs / shallow.Lp, 'shallow sea is steeper');
});

test('fine field: surface Hs matches the coarse estimate and the surface is bounded', () => {
  const erie = presetById('erie');
  const M = 600;
  const depths = sampleDepths(erie, M);
  const dx = erie.length / (M - 1);
  const wf = new WaveField();
  wf.computeCoarse(depths, dx, null, windStressFactor(15), +1, Infinity, { shoaling: true, breaking: true });
  const depthAt = depthFunction(erie);
  const x0 = 300e3, x1 = 300.5e3;
  const fine = wf.computeFine(x0, x1, 500, depthAt, null, { shoaling: true, breaking: true, nonlinear: true });
  // time-average the surface variance over many samples → Hs = 4 sqrt(m0)
  let sum = 0, cnt = 0;
  for (let t = 0; t < 600; t += 0.5) {
    const eta = wf.surface(t);
    for (let p = 0; p < eta.length; p += 10) { sum += eta[p] * eta[p]; cnt++; }
  }
  const hsMeasured = 4 * Math.sqrt(sum / cnt);
  const hsModel = fine.Hs[250];
  assert.ok(Math.abs(hsMeasured - hsModel) / hsModel < 0.2, `Hs measured ${hsMeasured} vs model ${hsModel}`);
});

test('slow model: converges to the equilibrium set-up and seiches at the long-wave period', () => {
  const erie = presetById('erie');
  const N = 200;
  const dx = erie.length / N;
  const depths = new Float64Array(N);
  const dfn = depthFunction(erie);
  for (let j = 0; j < N; j++) depths[j] = dfn((j + 0.5) * dx);
  const model = new SlowModel(depths, dx);
  model.reset();
  const tau = windStress(15);
  // blow for 96 hours
  for (let t = 0; t < 96 * 3600; t += 600) model.step(600, tau);
  const tilt = model.eta[N - 1] - model.eta[0];
  const eq = equilibriumSetup(model.d, dx, tau);
  const tiltEq = eq[N - 1] - eq[0];
  assert.ok(Math.abs(tilt - tiltEq) / tiltEq < 0.15, `tilt ${tilt} vs equilibrium ${tiltEq}`);
  assert.ok(tilt > 0.5 && tilt < 2.5, `Erie set-up ${tilt} m`);
  // stop the wind; measure the seiche period from upward zero crossings at the east end
  // (local maxima are unreliable: higher modes put small ripples on the fundamental)
  const crossings = [];
  let prev = model.eta[N - 1];
  for (let t = 0; t < 120 * 3600; t += 600) {
    model.step(600, 0);
    const e = model.eta[N - 1];
    if (prev < 0 && e >= 0) crossings.push(model.time);
    prev = e;
  }
  assert.ok(crossings.length >= 3, `found ${crossings.length} crossings`);
  const period = (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1);
  const estimate = seichePeriod(model.d, dx);
  assert.ok(Math.abs(period - estimate) / estimate < 0.15, `seiche ${period / 3600} h vs estimate ${estimate / 3600} h`);
  assert.ok(period / 3600 > 10 && period / 3600 < 18, `Erie seiche ${period / 3600} h`);
});

test('slow model with Erie widths: the narrow Buffalo end gains relative to a uniform channel', async () => {
  const { widthFunction } = await import('../src/bathymetry.js');
  const erie = presetById('erie');
  const N = 200;
  const dx = erie.length / N;
  const dfn = depthFunction(erie);
  const wfn = widthFunction(erie);
  const depths = new Float64Array(N);
  const widths = new Float64Array(N);
  for (let j = 0; j < N; j++) { depths[j] = dfn((j + 0.5) * dx); widths[j] = wfn((j + 0.5) * dx); }
  const tau = windStress(15);
  const run = (w) => {
    const model = new SlowModel(depths, dx, w);
    model.reset();
    for (let t = 0; t < 200 * 3600; t += 600) model.step(600, tau);
    // the impulsively started seiche rings for days: average over two periods
    let east = 0, west = 0, cnt = 0;
    const eta = new Float64Array(N);
    for (let t = 0; t < 32 * 3600; t += 600) {
      model.step(600, tau);
      east += model.eta[N - 1]; west += model.eta[0]; cnt++;
      for (let j = 0; j < N; j++) eta[j] += model.eta[j];
    }
    for (let j = 0; j < N; j++) eta[j] /= cnt;
    return { east: east / cnt, west: west / cnt, eta };
  };
  const uniform = run(null);
  const shaped = run(widths);
  assert.ok(shaped.east > 0 && shaped.west < 0, `east ${shaped.east} west ${shaped.west}`);
  // the wide, shallow western basin holds most of the volume, so the zero of the
  // tilted surface moves west and Buffalo's rise grows at Toledo's expense
  assert.ok(shaped.east > uniform.east, `Buffalo rise ${shaped.east} vs uniform ${uniform.east}`);
  assert.ok(-shaped.west < -uniform.west, `Toledo fall ${-shaped.west} vs uniform ${-uniform.west}`);
  // end-to-end tilt is set by depth alone and is unchanged by the width
  const tiltU = uniform.east - uniform.west;
  const tiltS = shaped.east - shaped.west;
  assert.ok(Math.abs(tiltS - tiltU) / tiltU < 0.05, `tilt ${tiltS} vs ${tiltU}`);
  // volume is conserved: Σ b η dx ≈ 0
  let vol = 0;
  for (let j = 0; j < N; j++) vol += widths[j] * shaped.eta[j];
  assert.ok(Math.abs(vol) < 1e-6 * widths.reduce((a, b) => a + b, 0), 'volume conserved');
});
