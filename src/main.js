// main.js — state, clocks, controls, readouts and the render loop.

import {
  G, RHO_WATER, TWO_PI, wavenumber, wavelength, deepWavelength, groupFactor, depthClass, ursell,
  windStressFactor, dragCoefficient, windStress, beaufort, windLabel, windLabelShort, smbShallow, smbDeep,
  durationLimitedFetch, fetchGrowthTime, jonswapShape, equilibriumSetup, seichePeriod, GAMMA_SIG,
  SMB_FITTED_MAX_WIND, CD_MAX,
} from './physics.js';
import { PRESETS, presetById, depthFunction, widthFunction, sampleDepths, cloneProfile, withLength } from './bathymetry.js';
import { WaveField } from './wavefield.js';
import { SlowModel } from './slowmodel.js';
import { VESSELS, vesselById, vesselPose, vesselShapes } from './vessel.js';
import { readPalette, fitCanvas, drawOverview, drawZoom, drawSpectrum, drawGauges, overviewGeometry, zoomGeometry, fmt } from './render.js';

const M_COARSE = 600;
const WIND_MAX = 231.5;        // m/s = 450 kn: far beyond any real wind, kept for play
const WIDTH_MIN = 100;         // m
const WIDTH_MAX = 10000;       // m
// the wind slider is quadratic so the everyday range keeps fine resolution
const windFromSlider = (v) => Math.round(WIND_MAX * (v / 100) * (v / 100) * 10) / 10;
const sliderFromWind = (U) => 100 * Math.sqrt(U / WIND_MAX);
const widthFromSlider = (v) => Math.round((WIDTH_MIN * Math.pow(WIDTH_MAX / WIDTH_MIN, v / 100)) / 10) * 10;
const sliderFromWidth = (w) => (100 * Math.log(w / WIDTH_MIN)) / Math.log(WIDTH_MAX / WIDTH_MIN);
const N_SLOW = 200;
const TAU_UP = 0.5 * 3600;     // lake seconds: how fast the sea responds to a rising wind
const TAU_DOWN = 2 * 3600;     // lake seconds: how fast a sea decays when the wind drops
const $ = (id) => document.getElementById(id);

const state = {
  presetId: 'erie',
  base: null,
  profile: null,
  edited: false,
  editing: false,
  U10: 15,
  dir: 1,
  lakeRate: 600,
  waveSpeed: 1,
  win: { center: 0, width: 400 },
  exag: 2,
  opts: { shoaling: true, breaking: true, nonlinear: true, setup: true, orbits: true },
  sea: { UAeff: 0, Fstate: 0 },
  waveTime: 0,
  lakeTime: 0,
  hoverX: null,
  compareId: '',
  vesselId: 'tanker',
};
let vesselShapeCache = { id: '', shapes: null };
let compare = null;   // { profile, name, depths, dx, wf, dir, coarse, hsMax, tiltEq, seicheT }

let L = 0;
let depthAt = null;
let widthAt = null;
let coarseDepths = null;
let coarseDx = 0;
let maxDepth = 1;
let slow = null;
let slowDx = 0;
const wf = new WaveField();
let pal = readPalette();
let dirty = true;
let particles = [];
let particleKey = '';
let frameNo = 0;
const setupScratch = new Float64Array(M_COARSE);

const ov = $('overview');
const zc = $('zoom');
const sc = $('spectrum');
const gc = $('gauges');
const tip = $('tip');

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

function windTau() { return windStress(state.U10) * state.dir; }

function buildModels({ keepWater }) {
  const prof = state.profile;
  L = prof.length;
  depthAt = depthFunction(prof);
  widthAt = widthFunction(prof);
  coarseDepths = sampleDepths(prof, M_COARSE);
  coarseDx = L / (M_COARSE - 1);
  maxDepth = 1;
  for (const d of coarseDepths) maxDepth = Math.max(maxDepth, d);
  slowDx = L / N_SLOW;
  const depths = new Float64Array(N_SLOW);
  const widths = new Float64Array(N_SLOW);
  for (let j = 0; j < N_SLOW; j++) {
    depths[j] = depthAt((j + 0.5) * slowDx);
    widths[j] = widthAt((j + 0.5) * slowDx);
  }
  if (slow && keepWater && slow.n === N_SLOW && Math.abs(slow.dx - slowDx) < 1e-9) {
    slow.setDepths(depths, widths);
  } else {
    slow = new SlowModel(depths, slowDx, widths);
    slow.reset(state.opts.setup ? equilibriumSetup(slow.d, slowDx, windTau(), slow.b) : null);
  }
  state.win.center = Math.max(state.win.width / 2, Math.min(L - state.win.width / 2, state.win.center));
  dirty = true;
  jump = true;
}

function buildCompare() {
  if (!state.compareId) { compare = null; return; }
  const profile = presetById(state.compareId);
  const M = 300;
  const depths = sampleDepths(profile, M);
  const dx = profile.length / (M - 1);
  const N = 120;
  const dfn = depthFunction(profile);
  const wfn = widthFunction(profile);
  const d = new Float64Array(N);
  const b = new Float64Array(N);
  for (let j = 0; j < N; j++) { d[j] = Math.max(0.5, dfn((j + 0.5) * profile.length / N)); b[j] = wfn((j + 0.5) * profile.length / N); }
  compare = { profile, name: profile.name.split(',')[0], depths, dx, wf: new WaveField(), d, b, slowDx: profile.length / N, coarse: null, hsMax: 0, tiltEq: 0, seicheT: seichePeriod(d, profile.length / N) };
}

function updateCompare() {
  if (!compare) return;
  const { sea, dir, opts } = state;
  compare.dir = dir;
  compare.coarse = compare.wf.computeCoarse(compare.depths, compare.dx, null, sea.UAeff, dir, sea.Fstate, opts);
  let hsMax = 0;
  for (let j = 0; j < compare.coarse.M; j++) hsMax = Math.max(hsMax, compare.coarse.Hs[j]);
  compare.hsMax = hsMax;
  const eq = equilibriumSetup(compare.d, compare.slowDx, windTau(), compare.b);
  compare.tiltEq = (eq[eq.length - 1] - eq[0]) * dir;
}

/** Hs and Tp on the comparison lake at the same distance from its upwind shore. */
function compareAtFetch(fetch) {
  const c = compare.coarse;
  const xm = compare.dir > 0 ? fetch : compare.profile.length - fetch;
  if (xm < 0 || xm > compare.profile.length) return null;
  const j = Math.max(0, Math.min(c.M - 1, Math.round(xm / compare.dx)));
  return { Hs: c.Hs[j], Tp: c.Tp[j], d: c.d[j] };
}

function steadySea() {
  state.sea.UAeff = windStressFactor(state.U10);
  state.sea.Fstate = L;
}

function advanceSea(dt) {
  const sea = state.sea;
  const UAt = windStressFactor(state.U10);
  const tau = UAt > sea.UAeff ? TAU_UP : TAU_DOWN;
  sea.UAeff += (UAt - sea.UAeff) * (1 - Math.exp(-dt / tau));
  if (sea.UAeff > 0.3) {
    const tEq = fetchGrowthTime(sea.UAeff, sea.Fstate);
    sea.Fstate = Math.min(L, durationLimitedFetch(sea.UAeff, tEq + dt));
  } else {
    sea.Fstate *= Math.exp(-dt / TAU_DOWN);
  }
}

function levelAt(x) { return state.opts.setup ? slow.levelAt(x) : 0; }

// Rebuild scheduling.  The field is rebuilt in two halves on separate frames
// (coarse march, then the window) so no single frame carries both, and the
// window's amplitude tables fade from the previous rebuild to the new one
// over the rebuild interval, so the sea evolves continuously between them.
const REBUILD_S = 0.4;
let jump = true;            // next fine rebuild replaces the tables outright (bed/direction/option change)
let coarseAt = -1;
let fineAt = -1;
let fineInterval = REBUILD_S;
let pendingFine = false;

function rebuildCoarse(nowS) {
  const { opts, sea, dir } = state;
  let setup = null;
  if (opts.setup) {
    for (let j = 0; j < M_COARSE; j++) setupScratch[j] = slow.levelAt(j * coarseDx);
    setup = setupScratch;
  }
  wf.computeCoarse(coarseDepths, coarseDx, setup, sea.UAeff, dir, sea.Fstate, opts);
  updateCompare();
  coarseAt = nowS;
  pendingFine = true;
}

function rebuildFine(nowS, blend) {
  const { opts, win } = state;
  const n = Math.max(200, Math.min(2400, Math.round(zc.clientWidth || 800)));
  const x0 = win.center - win.width / 2;
  const x1 = win.center + win.width / 2;
  wf.computeFine(x0, x1, n, depthAt, opts.setup ? (x) => slow.levelAt(x) : null, opts, blend);
  fineInterval = fineAt < 0 ? REBUILD_S : Math.max(0.05, Math.min(1, nowS - fineAt));
  fineAt = nowS;
  pendingFine = false;
  const key = `${x0}|${x1}|${n}|${state.exag}|${Math.round(hsMaxOf(wf.fine) / 2)}|${Math.round(wf.fine.eta0[Math.floor(n / 2)] / 2)}`;
  if (key !== particleKey) { particleKey = key; buildParticles(); }
}

/** Full rebuild right now (user action). */
function recompute(nowS = performance.now() / 1000) {
  rebuildCoarse(nowS);
  rebuildFine(nowS, !jump);
  jump = false;
}

function blendWeight(nowS) {
  return fineAt < 0 ? 1 : Math.min(1, (nowS - fineAt) / fineInterval);
}

// ---------------------------------------------------------------------------
// Particles (orbital motion)
// ---------------------------------------------------------------------------

const Z_LEVELS = [-0.75, -2, -4, -7, -11, -16, -23, -32, -45, -65, -90, -130, -180, -260];

function hsMaxOf(f) {
  let m = 0;
  for (let p = 0; p < f.n; p++) if (f.Hs[p] > m) m = f.Hs[p];
  return m;
}

function buildParticles() {
  particles = [];
  const f = wf.fine;
  if (!f) return;
  const H = zc.clientHeight || 300;
  const g = zoomGeometry(zc.clientWidth || 800, H, state.win.width, state.exag, hsMaxOf(f), f.eta0[Math.floor(f.n / 2)]);
  for (const frac of [0.1, 0.3, 0.5, 0.7, 0.9]) {
    const p = Math.round(frac * (f.n - 1));
    const d = f.d[p];
    if (d < 1) continue;
    for (const level of Z_LEVELS) {
      if (-level > 0.93 * d) break;
      const z = f.eta0[p] + level;              // below the local mean surface
      if (g.yOf(z) > H - 26) break;
      particles.push({ p, z, xRest: frac * state.win.width, trail: [], cur: [frac * state.win.width, z], lastPush: -1 });
    }
  }
}

function updateParticles(t, blendW) {
  for (const pt of particles) {
    const { dx, dz } = wf.particleDisplacement(pt.p, pt.z, t, blendW);
    pt.cur = [pt.xRest + dx, pt.z + dz];
    if (t - pt.lastPush >= 0.1 || t < pt.lastPush) {
      pt.lastPush = t;
      pt.trail.push(pt.cur);
      if (pt.trail.length > 70) pt.trail.shift();
    }
  }
}

// ---------------------------------------------------------------------------
// Controls
// ---------------------------------------------------------------------------

function setupControls() {
  const sel = $('preset');
  for (const p of PRESETS) {
    const o = document.createElement('option');
    o.value = p.id;
    o.textContent = p.name;
    sel.appendChild(o);
  }
  sel.value = state.presetId;
  sel.addEventListener('change', () => loadPreset(sel.value));

  const cmp = $('compare');
  for (const p of PRESETS) {
    const o = document.createElement('option');
    o.value = p.id;
    o.textContent = p.name;
    cmp.appendChild(o);
  }
  cmp.addEventListener('change', () => { state.compareId = cmp.value; buildCompare(); dirty = true; });

  const len = $('length');
  len.addEventListener('input', () => {
    const km = +len.value;
    $('lengthVal').textContent = `${km} km`;
    const frac = state.win.center / L;
    state.base = withLength(presetById(state.presetId), km * 1e3);
    state.profile = state.edited ? withLength(state.profile, km * 1e3) : cloneProfile(state.base);
    buildModels({ keepWater: false });
    state.win.center = frac * L;
    steadySea();
    refreshLabels();
  });

  $('editBed').addEventListener('click', () => {
    state.editing = !state.editing;
    $('editBed').classList.toggle('on', state.editing);
    $('editBed').setAttribute('aria-pressed', String(state.editing));
  });
  $('resetBed').addEventListener('click', () => {
    state.profile = cloneProfile(state.base);
    state.edited = false;
    buildModels({ keepWater: true });
    refreshLabels();
  });

  const wind = $('wind');
  wind.value = String(sliderFromWind(state.U10));
  wind.addEventListener('input', () => {
    state.U10 = windFromSlider(+wind.value);
    if (state.lakeRate === 0) steadySea();
    refreshLabels();
    dirty = true;
  });
  for (const b of $('dirSeg').querySelectorAll('button')) {
    b.addEventListener('click', () => {
      state.dir = +b.dataset.dir;
      for (const o of $('dirSeg').querySelectorAll('button')) o.classList.toggle('on', o === b);
      dirty = true;
      jump = true;
      refreshLabels();
    });
  }
  $('calmSea').addEventListener('click', () => {
    state.sea.UAeff = 0;
    state.sea.Fstate = 0;
    if (state.lakeRate === 0) { state.lakeRate = 600; $('lakeRate').value = '600'; }
    dirty = true;
  });

  $('lakeRate').addEventListener('change', () => {
    state.lakeRate = +$('lakeRate').value;
    if (state.lakeRate === 0) steadySea();
    dirty = true;
  });
  $('waveSpeed').addEventListener('change', () => { state.waveSpeed = +$('waveSpeed').value; refreshLabels(); });
  $('resetLake').addEventListener('click', () => { slow.reset(null); dirty = true; });

  const pos = $('winPos');
  pos.addEventListener('input', () => { setWindowCenter((+pos.value / 1000) * L); });
  const wid = $('winWidth');
  wid.value = String(sliderFromWidth(state.win.width));
  wid.addEventListener('input', () => {
    state.win.width = widthFromSlider(+wid.value);
    setWindowCenter(state.win.center);
  });
  const vs = $('vessel');
  for (const v of VESSELS) {
    const o = document.createElement('option');
    o.value = v.id;
    o.textContent = v.name;
    vs.appendChild(o);
  }
  vs.value = state.vesselId;
  vs.addEventListener('change', () => { state.vesselId = vs.value; });

  $('exag').addEventListener('change', () => { state.exag = +$('exag').value; dirty = true; refreshLabels(); });
  $('fitWaves').addEventListener('click', fitWaves);
  $('jumpUp').addEventListener('click', () => setWindowCenter(state.dir > 0 ? 1500 : L - 1500));
  $('jumpMid').addEventListener('click', () => setWindowCenter(L / 2));
  $('jumpDown').addEventListener('click', () => setWindowCenter(state.dir > 0 ? L - 1500 : 1500));

  const bind = (id, key) => {
    const el = $(id);
    el.addEventListener('change', () => {
      state.opts[key] = el.checked;
      if (key === 'setup' && el.checked) slow.reset(equilibriumSetup(slow.d, slowDx, windTau(), slow.b));
      dirty = true;
      jump = true;
    });
  };
  bind('optShoaling', 'shoaling');
  bind('optBreaking', 'breaking');
  bind('optNonlinear', 'nonlinear');
  bind('optSetup', 'setup');
  bind('optOrbits', 'orbits');
}

/** Widen or narrow the window to four peak wavelengths at its centre. */
function fitWaves() {
  const s = windowCenterStats();
  if (!s.wet || !(s.Lp > 0)) return;
  setWindowWidth(4 * s.Lp);
}

function setWindowWidth(w) {
  state.win.width = Math.round(Math.max(WIDTH_MIN, Math.min(WIDTH_MAX, w)) / 10) * 10;
  $('winWidth').value = String(sliderFromWidth(state.win.width));
  setWindowCenter(state.win.center);
}

function setWindowCenter(x) {
  const half = state.win.width / 2;
  state.win.center = Math.max(half, Math.min(L - half, x));
  $('winPos').value = String(Math.round((state.win.center / L) * 1000));
  dirty = true;
  refreshLabels();
}

function loadPreset(id) {
  state.presetId = id;
  const p = presetById(id);
  state.base = cloneProfile(p);
  state.profile = cloneProfile(p);
  state.edited = false;
  $('length').value = String(Math.round(p.length / 1e3));
  $('lengthVal').textContent = `${Math.round(p.length / 1e3)} km`;
  buildModels({ keepWater: false });
  state.win.center = state.dir > 0 ? L - 8000 : 8000;
  steadySea();
  setWindowCenter(state.win.center);
  refreshLabels();
}

function refreshLabels() {
  const p = state.profile;
  const U = state.U10;
  $('windVal').textContent = `${fmt(U, 1)} m/s · ${fmt(U * 1.944, 0)} kn · ${windLabelShort(U)}`;
  $('dirWest').textContent = `◀ ${p.ends[0]}`;
  $('dirEast').textContent = `${p.ends[1]} ▶`;
  $('profileNote').textContent = `${p.stats}${state.edited ? ' · edited' : ''}. ${p.note}`;
  $('winPosVal').textContent = `${fmt(state.win.center / 1e3, 1)} km from ${p.ends[0]}`;
  $('winWidthVal').textContent = `${state.win.width} m`;
  $('zoomHint').textContent = `${state.win.width} m wide, vertical ×${state.exag} (ship included), waves at ${state.waveSpeed === 0 ? 'rest' : state.waveSpeed + '× real time'}`;
  $('legWest').textContent = p.ends[0];
  $('legEast').textContent = p.ends[1];
  $('jumpUp').textContent = `${state.dir > 0 ? p.ends[0] : p.ends[1]} (upwind)`;
  $('jumpDown').textContent = `${state.dir > 0 ? p.ends[1] : p.ends[0]} (downwind)`;
}

// ---------------------------------------------------------------------------
// Expand a panel: real fullscreen where the browser allows it, otherwise an
// in-page overlay that fills the viewport.  Escape or the button closes it.
// ---------------------------------------------------------------------------

let maxPanel = null;

function setupExpand() {
  for (const btn of document.querySelectorAll('[data-expand]')) {
    const panel = btn.closest('.panel');
    btn.addEventListener('click', () => (maxPanel === panel ? unmaximize() : maximize(panel)));
  }
  document.addEventListener('fullscreenchange', () => {
    if (!document.fullscreenElement && maxPanel && maxPanel.dataset.fs === '1') unmaximize();
  });
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && maxPanel) { unmaximize(); return; }
    const t = e.target;
    const tag = t && t.tagName;
    const inSlider = tag === 'INPUT' && t.type === 'range';
    const inField = tag === 'TEXTAREA' || tag === 'SELECT' || (tag === 'INPUT' && !inSlider && t.type !== 'checkbox' && t.type !== 'button');
    if (inField) return;
    if ((e.key === 'f' || e.key === 'F') && !e.ctrlKey && !e.metaKey && !e.altKey) { fitWaves(); return; }
    if (inSlider) return;                       // arrows belong to the focused slider
    const big = e.shiftKey;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      const step = big ? 1000 : state.win.width * 0.25;
      setWindowCenter(state.win.center + (e.key === 'ArrowRight' ? step : -step));
      e.preventDefault();
    } else if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
      const factor = big ? 2 : 1.25;
      setWindowWidth(e.key === 'ArrowUp' ? state.win.width / factor : state.win.width * factor);
      e.preventDefault();
    }
  });
}

function setExpandLabels() {
  for (const btn of document.querySelectorAll('[data-expand]')) {
    const on = btn.closest('.panel') === maxPanel;
    btn.textContent = on ? 'Close (Esc)' : 'Expand';
    btn.setAttribute('aria-pressed', String(on));
  }
}

function maximize(panel) {
  if (maxPanel) unmaximize();
  maxPanel = panel;
  panel.classList.add('maximized');
  document.body.classList.add('has-max');
  setExpandLabels();
  if (panel.requestFullscreen) {
    panel.requestFullscreen().then(() => { panel.dataset.fs = '1'; }).catch(() => { /* overlay fallback is already in place */ });
  }
  dirty = true;
}

function unmaximize() {
  const panel = maxPanel;
  if (!panel) return;
  maxPanel = null;
  panel.classList.remove('maximized');
  document.body.classList.remove('has-max');
  delete panel.dataset.fs;
  if (document.fullscreenElement === panel && document.exitFullscreen) document.exitFullscreen().catch(() => {});
  setExpandLabels();
  dirty = true;
}

// ---------------------------------------------------------------------------
// Overview interaction: drag the window, drag bed points, hover for values
// ---------------------------------------------------------------------------

function setupOverviewInteraction() {
  let drag = null;   // { kind: 'window' | 'point', index }
  const geom = () => overviewGeometry(ov.clientWidth, ov.clientHeight, L, maxDepth);
  const rel = (e) => { const r = ov.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };

  ov.addEventListener('pointerdown', (e) => {
    const { x, y } = rel(e);
    const g = geom();
    if (state.editing) {
      const pts = state.profile.points;
      let best = -1;
      let bestD = 12;
      for (let i = 1; i < pts.length - 1; i++) {
        const px = g.xOf(pts[i][0] * 1e3);
        const py = g.yOfDepth(pts[i][1]);
        const dd = Math.hypot(px - x, py - y);
        if (dd < bestD) { bestD = dd; best = i; }
      }
      if (best >= 0) { drag = { kind: 'point', index: best }; ov.setPointerCapture(e.pointerId); return; }
    }
    drag = { kind: 'window' };
    ov.setPointerCapture(e.pointerId);
    setWindowCenter(g.xmOf(x));
    tip.hidden = true;
  });
  ov.addEventListener('pointermove', (e) => {
    const { x, y } = rel(e);
    const g = geom();
    if (drag && drag.kind === 'window') { setWindowCenter(g.xmOf(x)); return; }
    if (drag && drag.kind === 'point') {
      const d = Math.max(0, Math.min(600, (y - g.yWater) / g.pxPerM));
      state.profile.points[drag.index][1] = Math.round(d * 10) / 10;
      state.edited = true;
      buildModels({ keepWater: true });
      refreshLabels();
      return;
    }
    const xm = g.xmOf(x);
    if (xm < 0 || xm > L || y < g.chartTop || y > g.secBottom) { state.hoverX = null; tip.hidden = true; return; }
    state.hoverX = xm;
    showTip(e.clientX, e.clientY, xm);
  });
  const end = () => { drag = null; };
  ov.addEventListener('pointerup', end);
  ov.addEventListener('pointercancel', end);
  ov.addEventListener('pointerleave', () => { state.hoverX = null; tip.hidden = true; });
}

function showTip(cx, cy, xm) {
  const c = wf.coarse;
  if (!c) return;
  const j = Math.max(0, Math.min(c.M - 1, Math.round(xm / coarseDx)));
  const d = c.d[j];
  const Tp = c.Tp[j];
  const lines = [`${fmt(xm / 1e3, 1)} km · depth ${fmt(d, 1)} m`];
  if (c.Hs[j] > 0.005 && d > 0.3) {
    const k = wavenumber(TWO_PI / Tp, d);
    lines.push(`Hs ${fmt(c.Hs[j], 2)} m · Tp ${fmt(Tp, 1)} s · L ${fmt(TWO_PI / k, 0)} m`);
    lines.push(`d/L ${fmt(d / (TWO_PI / k), 2)} (${depthClass(k * d)})${c.breaking[j] ? ' · breaking' : ''}`);
  }
  if (state.opts.setup) lines.push(`level ${levelAt(xm) >= 0 ? '+' : ''}${fmt(levelAt(xm), 2)} m`);
  tip.innerHTML = lines.map((s) => `<div>${s}</div>`).join('');
  tip.hidden = false;
  const r = tip.getBoundingClientRect();
  const left = Math.min(window.innerWidth - r.width - 8, cx + 14);
  tip.style.left = `${Math.max(4, left)}px`;
  tip.style.top = `${Math.max(4, cy - r.height - 12)}px`;
}

// ---------------------------------------------------------------------------
// Readouts
// ---------------------------------------------------------------------------

const READOUT_SPEC = [
  { title: 'Wind', rows: [['u10', 'Speed at 10 m'], ['ua', 'Stress factor U_A'], ['cd', 'Drag coefficient C_D'], ['tau', 'Wind stress τ']] },
  { title: 'Sea growth', rows: [['fetch', 'Fetch to downwind shore'], ['age', 'Fetch the wind has filled'], ['hsEnd', 'Hs at downwind shore'], ['tpEnd', 'Tp at downwind shore'], ['deepEq', 'Same wind over deep water']] },
  { title: 'At the window', rows: [['depth', 'Depth'], ['hs', 'Hs'], ['tp', 'Tp'], ['lp', 'Peak wavelength L'], ['regime', 'd / L'], ['c', 'Phase speed C'], ['cg', 'Group speed Cg'], ['ks', 'Shoaling Ks'], ['steep', 'Steepness Hs / L'], ['hd', 'Hs / d'], ['ur', 'Ursell number']] },
  { title: 'Whole basin', rows: [['meanD', 'Mean depth along section'], ['maxD', 'Max depth'], ['lakeStats', 'Lake, area-wide'], ['seiche', 'Seiche period T₁'], ['tilt', 'Tilt now, downwind − upwind'], ['tiltEq', 'Equilibrium tilt'], ['clock', 'Lake clock']] },
  { title: 'Compared lake, same wind', rows: [['cmpName', 'Lake'], ['cmpHere', 'Hs · Tp at the same fetch'], ['cmpDepth', 'Depth there'], ['cmpEnd', 'Max Hs along its fetch'], ['cmpTilt', 'Equilibrium tilt'], ['cmpSeiche', 'Seiche period T₁']] },
];
const readoutEls = {};

function buildReadouts() {
  const root = $('readouts');
  for (const card of READOUT_SPEC) {
    const el = document.createElement('div');
    el.className = 'card';
    const h = document.createElement('h3');
    h.textContent = card.title;
    el.appendChild(h);
    const dl = document.createElement('dl');
    for (const [id, label] of card.rows) {
      const dt = document.createElement('dt');
      dt.textContent = label;
      const dd = document.createElement('dd');
      dd.id = `ro-${id}`;
      dd.textContent = '–';
      dl.appendChild(dt);
      dl.appendChild(dd);
      readoutEls[id] = dd;
    }
    el.appendChild(dl);
    root.appendChild(el);
  }
}

function clockText(seconds) {
  const h = Math.floor(seconds / 3600);
  const m = Math.floor((seconds % 3600) / 60);
  return `${h} h ${String(m).padStart(2, '0')} min`;
}

function windowCenterStats() {
  const f = wf.fine;
  const p = Math.floor(f.n / 2);
  const d = f.d[p];
  const Hs = f.Hs[p];
  const Tp = f.Tp[p];
  const wet = d > 0.3 && Hs > 1e-4;
  const k = wet ? wavenumber(TWO_PI / Tp, d) : 0;
  const Lp = wet ? TWO_PI / k : 0;
  const C = wet ? Lp / Tp : 0;
  const Cg = wet ? groupFactor(k * d) * C : 0;
  const dEff = f.dEff[p];
  const kG = wet ? wavenumber(TWO_PI / Tp, dEff) : 0;
  const CgGen = wet ? groupFactor(kG * dEff) * (TWO_PI / kG) / Tp : 0;
  const Ks = wet ? Math.sqrt(CgGen / Cg) : 0;
  return { p, d, Hs, Tp, k, Lp, C, Cg, Ks, dEff, wet, breaking: !!f.breaking[p], L0: wet ? deepWavelength(Tp) : 0, C0: wet ? deepWavelength(Tp) / Tp : 0 };
}

function updateReadouts() {
  const U = state.U10;
  const { sea, profile } = state;
  const c = wf.coarse;
  const s = windowCenterStats();
  const set = (id, text, hot = false) => { const el = readoutEls[id]; el.textContent = text; el.classList.toggle('hot', hot); };
  set('u10', `${fmt(U, 1)} m/s · ${fmt(U * 1.944, 0)} kn · ${fmt(U * 2.237, 0)} mph`, U > SMB_FITTED_MAX_WIND);
  set('ua', `${fmt(sea.UAeff, 1)} m/s (target ${fmt(windStressFactor(U), 1)})`);
  set('cd', `${(dragCoefficient(U) * 1e3).toFixed(2)} × 10⁻³`);
  set('tau', `${fmt(windStress(U), 3)} Pa`);
  set('fetch', `${fmt(L / 1e3, 0)} km`);
  const Feff = Math.min(L, sea.Fstate);
  set('age', sea.UAeff < 0.3 ? 'no wind' : sea.Fstate >= L - 1 ? 'all of it (fetch-limited)' : `${fmt(Feff / 1e3, 0)} km (still growing)`);
  const endJ = state.dir > 0 ? c.M - 1 : 0;
  let hsMax = 0;
  let jMax = endJ;
  for (let j = 0; j < c.M; j++) if (c.HsGen[j] > hsMax) { hsMax = c.HsGen[j]; jMax = j; }
  set('hsEnd', `${fmt(c.HsGen[jMax], 2)} m`);
  set('tpEnd', `${fmt(c.Tp[jMax], 1)} s`);
  const deep = smbDeep(sea.UAeff, Feff);
  set('deepEq', `${fmt(deep.Hs, 2)} m, ${fmt(deep.Tp, 1)} s`);
  set('depth', `${fmt(s.d, 1)} m`);
  set('hs', `${fmt(s.Hs, 2)} m`, s.breaking);
  set('tp', `${fmt(s.Tp, 1)} s`);
  set('lp', s.wet ? `${fmt(s.Lp, 0)} m (deep: ${fmt(s.L0, 0)})` : '–');
  set('regime', s.wet ? `${fmt(s.d / s.Lp, 2)} · ${depthClass(s.k * s.d)}` : '–');
  set('c', s.wet ? `${fmt(s.C, 1)} m/s (deep: ${fmt(s.C0, 1)})` : '–');
  set('cg', s.wet ? `${fmt(s.Cg, 1)} m/s` : '–');
  set('ks', s.wet ? `${fmt(s.Ks, 2)} from ${fmt(s.dEff, 0)} m` : '–');
  set('steep', s.wet ? fmt(s.Hs / s.Lp, 3) : '–');
  set('hd', s.wet ? `${fmt(s.Hs / s.d, 2)}${s.breaking ? ' · breaking' : ` (limit ${GAMMA_SIG})`}` : '–', s.breaking);
  set('ur', s.wet ? fmt(ursell(s.Hs, s.Lp, s.d), 1) : '–');
  let meanD = 0;
  for (let j = 0; j < c.M; j++) meanD += coarseDepths[j];
  meanD /= c.M;
  set('meanD', `${fmt(meanD, 0)} m`);
  set('maxD', `${fmt(maxDepth, 0)} m`);
  set('lakeStats', profile.stats);
  set('seiche', `${fmt(seichePeriod(slow.d, slowDx) / 3600, 1)} h`);
  const tiltNow = (slow.eta[slow.n - 1] - slow.eta[0]) * state.dir;
  set('tilt', state.opts.setup ? `${fmt(tiltNow, 2)} m` : 'off');
  const eq = equilibriumSetup(slow.d, slowDx, windTau(), slow.b);
  set('tiltEq', `${fmt((eq[slow.n - 1] - eq[0]) * state.dir, 2)} m`);
  set('clock', clockText(state.lakeTime));
  $('clock').textContent = `lake clock ${clockText(state.lakeTime)}`;
  if (compare && compare.coarse) {
    const fetchHere = state.dir > 0 ? state.win.center : L - state.win.center;
    const at = compareAtFetch(fetchHere);
    set('cmpName', compare.name);
    set('cmpHere', at ? `${fmt(at.Hs, 2)} m · ${fmt(at.Tp, 1)} s` : 'beyond its length');
    set('cmpDepth', at ? `${fmt(at.d, 0)} m` : '–');
    set('cmpEnd', `${fmt(compare.hsMax, 2)} m`);
    set('cmpTilt', `${fmt(compare.tiltEq, 2)} m`);
    set('cmpSeiche', `${fmt(compare.seicheT / 3600, 1)} h`);
  } else {
    for (const id of ['cmpName', 'cmpHere', 'cmpDepth', 'cmpEnd', 'cmpTilt', 'cmpSeiche']) set(id, '–');
    set('cmpName', 'choose one under Lake');
  }
}

// ---------------------------------------------------------------------------
// "Why" text
// ---------------------------------------------------------------------------

function currentVessel(eta) {
  const def = vesselById(state.vesselId);
  if (!def.LOA) return null;
  if (vesselShapeCache.id !== def.id) vesselShapeCache = { id: def.id, shapes: vesselShapes(def) };
  const f = wf.fine;
  const pc = Math.floor(f.n / 2);
  if (f.d[pc] < 0.3) return null;
  const pose = vesselPose(def, f, eta, pc, -state.dir);   // bow into the wind
  return { def, pose, shapes: vesselShapeCache.shapes };
}

function updateWhy() {
  const { sea, profile, dir, opts } = state;
  const U = state.U10;
  const c = wf.coarse;
  const up = dir > 0 ? profile.ends[0] : profile.ends[1];
  const down = dir > 0 ? profile.ends[1] : profile.ends[0];
  const paras = [];
  if (sea.UAeff < 0.5 || U < 0.5) {
    paras.push(`With no wind there is nothing to raise waves. Set a wind speed, or wait for the sea to build on the lake clock after a calm.`);
  } else {
    const Feff = Math.min(L, sea.Fstate);
    const endJ = dir > 0 ? c.M - 1 : 0;
    let jBest = endJ;
    for (let j = 0; j < c.M; j++) if (c.HsGen[j] > c.HsGen[jBest]) jBest = j;
    const dEff = c.dEff[jBest];
    const here = smbShallow(sea.UAeff, Feff, dEff);
    const deep = smbDeep(sea.UAeff, Feff);
    const steepHere = here.Hs / wavelength(here.Tp, dEff);
    const steepDeep = deep.Hs / deepWavelength(deep.Tp);
    let p1 = `A ${fmt(U, 0)} m/s wind (${fmt(U * 1.944, 0)} kn, ${windLabel(U)}) blowing toward ${down} works on ${fmt(Feff / 1e3, 0)} km of water that averages ${fmt(dEff, 0)} m deep along the way. The growth curves give a sea of Hs ${fmt(here.Hs, 1)} m and Tp ${fmt(here.Tp, 1)} s at the ${down} end.`;
    if (sea.Fstate < L - 1) p1 += ` The wind has not blown long enough to use the whole ${fmt(L / 1e3, 0)} km fetch, so the sea is duration-limited and still growing on the lake clock.`;
    if (deep.Hs > here.Hs * 1.03) {
      p1 += ` Over deep water the same wind and fetch would raise ${fmt(deep.Hs, 1)} m waves of ${fmt(deep.Tp, 1)} s. The bed drains energy from the longest waves, so this sea is lower and shorter, and steeper: Hs/L is ${fmt(steepHere, 3)} here against ${fmt(steepDeep, 3)} in deep water. Short, steep waves are what make a shallow lake dangerous in a blow.`;
    } else {
      p1 += ` That is essentially the deep-water result: the water along the fetch is deep enough that the bed never limits growth.`;
    }
    if (U > SMB_FITTED_MAX_WIND) {
      p1 += ` Above about ${SMB_FITTED_MAX_WIND} m/s the growth curves are extrapolated far beyond the storms they were fitted to, and the drag coefficient is held at its hurricane ceiling of ${(CD_MAX * 1e3).toFixed(1)} × 10⁻³: treat these numbers as play, not prediction.`;
    }
    paras.push(p1);

    const s = windowCenterStats();
    if (s.wet) {
      const cls = depthClass(s.k * s.d);
      let p2 = `At the window the depth is ${fmt(s.d, 1)} m, which is ${fmt(s.d / s.Lp, 2)} wavelengths: ${cls} water.`;
      if (cls === 'deep') {
        p2 += ` The bed is more than half a wavelength down, so the ${fmt(s.Tp, 1)} s waves do not feel it. They travel at the deep-water speed of ${fmt(s.C0, 1)} m/s with a wavelength of ${fmt(s.L0, 0)} m, and the orbits are circles that fade out well above the bottom.`;
      } else {
        p2 += ` The ${fmt(s.Tp, 1)} s waves feel the bottom: compared with deep water they have slowed from ${fmt(s.C0, 1)} to ${fmt(s.C, 1)} m/s and shortened from ${fmt(s.L0, 0)} to ${fmt(s.Lp, 0)} m. The orbits flatten into ellipses${cls === 'shallow' ? ' that reach the bed, which is why a storm stirs up the mud' : ''}.`;
        if (opts.shoaling) p2 += ` Energy flux is conserved, so the height scales by Ks = ${fmt(s.Ks, 2)} relative to the ${fmt(s.dEff, 0)} m water it grew in.`;
      }
      const wavesInWindow = state.win.width / s.Lp;
      if (wavesInWindow < 2) {
        p2 += ` The window is only ${fmt(wavesInWindow, 1)} wavelengths wide, so instead of waves marching past you see the whole surface heave and tilt as each ${fmt(s.Lp, 0)} m crest passes; the short chop on top is the high-frequency part of the same wind-wave spectrum. Press F or widen the window to see several waves at once.`;
      }
      if (opts.breaking) {
        p2 += s.breaking
          ? ` Hs/d has reached ${fmt(s.Hs / s.d, 2)}, past the ${GAMMA_SIG} limit, so the waves are breaking and the surplus energy goes into foam.`
          : ` Hs/d is ${fmt(s.Hs / s.d, 2)}, below the ${GAMMA_SIG} breaking limit.`;
      }
      const def = vesselById(state.vesselId);
      if (def.LOA) {
        const ratio = def.LOA / s.Lp;
        const h = s.d - def.draft;
        let p2b = `The ${def.name.toLowerCase()} drawn in the window has a draft of ${def.draft} m, leaving ${fmt(h, 1)} m of water under the keel here${h < 0 ? ': it is aground' : ''}.`;
        if (ratio > 1.8) p2b += ` Its ${def.LOA} m hull spans ${fmt(ratio, 1)} wavelengths, so the crests and troughs average out along it and it rides nearly level.`;
        else if (ratio > 0.6) p2b += ` Its ${def.LOA} m hull is close to one wavelength (${fmt(ratio, 1)}): this is the sea that pitches a ship hardest and alternately hogs and sags the hull.`;
        else p2b += ` At ${fmt(ratio, 2)} wavelengths long it is small compared with the waves and simply rides over them.`;
        paras.push(p2 + ' ' + p2b);
      } else {
        paras.push(p2);
      }
    } else {
      paras.push(`The window is on dry ground or over water too shallow to carry waves. Move it along the lake.`);
    }
  }
  const tau = windStress(U);
  const eq = equilibriumSetup(slow.d, slowDx, tau * dir, slow.b);
  const tiltEq = (eq[slow.n - 1] - eq[0]) * dir;
  const T1 = seichePeriod(slow.d, slowDx) / 3600;
  let meanD = 0;
  for (let j = 0; j < slow.n; j++) meanD += slow.d[j];
  meanD /= slow.n;
  const compareDepth = meanD < 60 ? 150 : 15;
  const tiltCompare = (tau * L) / (RHO_WATER * G * compareDepth);
  let p3 = `The whole lake responds too. A stress of ${fmt(tau, 2)} Pa on the surface tilts it until the slope of the water balances the push: dη/dx = τ/(ρgh), so the shallower the water, the steeper the tilt.`;
  if (tau > 0) {
    p3 += ` Along this section the equilibrium set-up is ${fmt(tiltEq, 2)} m higher at ${down} than at ${up}` + (opts.setup ? `, and the model currently shows ${fmt((slow.eta[slow.n - 1] - slow.eta[0]) * dir, 2)} m.` : '.');
    p3 += ` A basin of the same length but ${compareDepth} m deep would tilt only about ${fmt(tiltCompare, 2)} m under the same wind.`;
  }
  if (compare && compare.coarse && tau > 0) {
    p3 += ` ${compare.name} under this wind: equilibrium tilt ${fmt(compare.tiltEq, 2)} m, seiche period ${fmt(compare.seicheT / 3600, 1)} h, biggest waves ${fmt(compare.hsMax, 1)} m.`;
  }
  p3 += ` If the wind changes, the tilted surface overshoots and rocks back as a seiche with a period near ${fmt(T1, 1)} h: the time for a long wave moving at √(gh) to cross the basin and return.`;
  paras.push(p3);
  $('why').innerHTML = paras.map((p) => `<p>${p}</p>`).join('');
}

// ---------------------------------------------------------------------------
// Charts
// ---------------------------------------------------------------------------

function renderSpectrum() {
  const { ctx, W, H } = fitCanvas(sc);
  const f = wf.fine;
  const p = Math.floor(f.n / 2);
  const N = wf.nComp;
  const sHere = new Float64Array(N);
  const sDeep = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const a = f.amp[i * f.n + p];
    sHere[i] = (a * a) / (2 * wf.df[i]);
  }
  const Feff = Math.min(L, state.sea.Fstate);
  const deep = smbDeep(state.sea.UAeff, Feff);
  let m0 = 0;
  for (let i = 0; i < N; i++) { sDeep[i] = deep.Tp > 0 ? jonswapShape(wf.f[i], 1 / deep.Tp) * wf.df[i] : 0; m0 += sDeep[i]; }
  const scale = m0 > 0 ? (deep.Hs * deep.Hs) / 16 / m0 : 0;
  for (let i = 0; i < N; i++) sDeep[i] = (sDeep[i] * scale) / wf.df[i];
  drawSpectrum(ctx, W, H, { f: wf.f, sHere, sDeep, fpHere: f.Tp[p] > 0 ? 1 / f.Tp[p] : 0, fpDeep: deep.Tp > 0 ? 1 / deep.Tp : 0, hsHere: f.Hs[p], hsDeep: deep.Hs }, pal);
}

function renderGauges() {
  const { ctx, W, H } = fitCanvas(gc);
  drawGauges(ctx, W, H, { samples: slow.samples, now: slow.time, hours: 48, ends: state.profile.ends, seicheT: seichePeriod(slow.d, slowDx) }, pal);
}

// ---------------------------------------------------------------------------
// Loop
// ---------------------------------------------------------------------------

let lastNow = performance.now();

function frame(now) {
  const nowS = now / 1000;
  const dtReal = Math.min(0.1, Math.max(0, (now - lastNow) / 1000));
  lastNow = now;
  frameNo++;
  if (state.waveSpeed > 0) state.waveTime += dtReal * state.waveSpeed;
  if (state.lakeRate > 0) {
    const dtL = dtReal * state.lakeRate;
    advanceSea(dtL);
    if (state.opts.setup) slow.step(dtL, windTau());
    state.lakeTime += dtL;
  }
  if (dirty && nowS - fineAt >= 0.08) {
    recompute(nowS);
    dirty = false;
  } else if (state.lakeRate > 0) {
    if (pendingFine) rebuildFine(nowS, true);
    else if (nowS - coarseAt >= REBUILD_S) rebuildCoarse(nowS);
  }
  const blendW = blendWeight(nowS);

  // overview
  {
    const { ctx, W, H } = fitCanvas(ov);
    drawOverview(ctx, W, H, {
      L, coarse: wf.coarse, depthAt, levelAt: state.opts.setup ? (x) => slow.levelAt(x) : null,
      win: state.win, dir: state.dir, U10: state.U10, ends: state.profile.ends,
      editing: state.editing, points: state.profile.points, maxDepth, hoverX: state.hoverX,
      compare: compare && compare.coarse ? compare : null,
    }, pal);
  }
  // zoom
  {
    const { ctx, W, H } = fitCanvas(zc);
    const eta = wf.surface(state.waveTime, blendW);
    if (state.opts.orbits) updateParticles(state.waveTime, blendW);
    const s = windowCenterStats();
    drawZoom(ctx, W, H, {
      fine: wf.fine, eta, win: state.win, exag: state.exag, particles, dir: state.dir, U10: state.U10,
      L, ends: state.profile.ends, opts: state.opts, center: state.win.center, Lp: s.Lp,
      vessel: currentVessel(eta), hsMax: hsMaxOf(wf.fine),
    }, pal);
  }
  if (frameNo % 6 === 0) { renderSpectrum(); renderGauges(); }
  if (frameNo % 12 === 0) updateReadouts();
  if (frameNo % 40 === 0) updateWhy();
  requestAnimationFrame(frame);
}

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------

function boot() {
  setupControls();
  buildReadouts();
  loadPreset(state.presetId);
  setupOverviewInteraction();
  setupExpand();
  const mq = window.matchMedia('(prefers-color-scheme: dark)');
  const refreshPal = () => { pal = readPalette(); };
  if (mq.addEventListener) mq.addEventListener('change', refreshPal);
  new MutationObserver(refreshPal).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  new ResizeObserver(() => { dirty = true; }).observe(zc);
  recompute();
  dirty = false;
  updateReadouts();
  updateWhy();
  requestAnimationFrame(frame);
}

boot();
