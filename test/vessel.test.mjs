import test from 'node:test';
import assert from 'node:assert/strict';
import { VESSELS, fitWaterline, vesselPose, vesselShapes } from '../src/vessel.js';

test('waterline fit recovers a tilted plane surface', () => {
  const n = 401;
  const x = Float64Array.from({ length: n }, (_, p) => p * 0.5);
  const eta = Float64Array.from(x, (xx) => 0.4 + 0.02 * (xx - 100));
  const fit = fitWaterline(eta, x, 200, 75);
  assert.ok(Math.abs(fit.b - 0.4) < 1e-9, `b ${fit.b}`);
  assert.ok(Math.abs(fit.m - 0.02) < 1e-9, `m ${fit.m}`);
  assert.equal(fit.samples, 301);
});

test('a long hull filters short waves, a short hull rides them', () => {
  const n = 2001;
  const x = Float64Array.from({ length: n }, (_, p) => p * 0.25);  // 500 m window
  const Lw = 60;                                                   // 60 m waves, 2 m amplitude
  const eta = Float64Array.from(x, (xx) => 2 * Math.cos((2 * Math.PI * xx) / Lw));
  const pc = 960;                                  // x = 240 m, exactly four wavelengths: a crest
  const long = fitWaterline(eta, x, pc, 75);      // 150 m hull: 2.5 wavelengths
  const short = fitWaterline(eta, x, pc, 4);      // 8 m boat
  assert.ok(Math.abs(long.b) < 0.4, `long hull heave ${long.b}`);
  assert.ok(Math.abs(short.b - 2) < 0.1, `boat rides the crest: ${short.b}`);
});

test('vessel pose: afloat when the water is deeper than the draft, aground otherwise', () => {
  const def = VESSELS.find((v) => v.id === 'tanker');
  const n = 201;
  const x = Float64Array.from({ length: n }, (_, p) => p * 2);
  const fine = { x, d: new Float64Array(n).fill(20), eta0: new Float64Array(n) };
  const eta = new Float64Array(n).fill(0.5);
  const afloat = vesselPose(def, fine, eta, 100, 1);
  assert.equal(afloat.aground, false);
  assert.ok(Math.abs(afloat.heave - 0.5) < 1e-9);
  fine.d.fill(6);
  const stuck = vesselPose(def, fine, eta, 100, 1);
  assert.equal(stuck.aground, true);
  assert.ok(Math.abs(stuck.heave - (8 - 6)) < 1e-9, `keel on the bed: heave ${stuck.heave}`);
});

test('every vessel has a hull polygon spanning its length and draft', () => {
  for (const v of VESSELS) {
    if (v.id === 'none') continue;
    const shapes = vesselShapes(v);
    const hull = shapes.find((s) => s.part === 'hull');
    const xs = hull.pts.map((p) => p[0]);
    const zs = hull.pts.map((p) => p[1]);
    assert.ok(Math.abs(Math.max(...xs) - Math.min(...xs) - v.LOA) < 1e-9, `${v.id} length`);
    assert.ok(Math.abs(Math.min(...zs) + v.draft) < 1e-9, `${v.id} draft`);
    const freeboard = v.hullDepth - v.draft;
    assert.ok(Math.max(...zs) >= freeboard - 1e-9 && Math.max(...zs) <= freeboard + 0.5, `${v.id} freeboard (sheer allowed)`);
  }
});
