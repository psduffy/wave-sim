// render.js — canvas drawing for the overview, the zoom window and the two charts.
// Colours come from the page's CSS tokens so the canvases follow the theme.

const TOKENS = ['fg', 'muted', 'line', 'accent', 'accent-soft', 'water-top', 'water-deep', 'water-fill', 'bed', 'bed-dark', 'sky', 'foam', 'breaking', 'wind', 'west', 'east', 'grid', 'panel', 'ship-hull', 'ship-edge', 'ship-boot', 'ship-house', 'ship-funnel', 'ship-fitting'];

export function readPalette() {
  const cs = getComputedStyle(document.documentElement);
  const pal = {};
  for (const t of TOKENS) pal[t.replace(/-([a-z])/g, (_, c) => c.toUpperCase())] = cs.getPropertyValue('--' + t).trim();
  pal.mono = cs.getPropertyValue('--font-mono').trim() || 'monospace';
  pal.sans = cs.getPropertyValue('--font-sans').trim() || 'sans-serif';
  return pal;
}

/** Size the backing store to the CSS box × devicePixelRatio; returns CSS-px dims. */
export function fitCanvas(canvas) {
  const dpr = Math.min(3, window.devicePixelRatio || 1);
  const W = Math.max(1, Math.round(canvas.clientWidth));
  const H = Math.max(1, Math.round(canvas.clientHeight));
  const bw = Math.round(W * dpr);
  const bh = Math.round(H * dpr);
  if (canvas.width !== bw || canvas.height !== bh) {
    canvas.width = bw;
    canvas.height = bh;
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W, H };
}

export function niceStep(span, targetCount) {
  const raw = span / Math.max(1, targetCount);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= raw) return m * p;
  return 10 * p;
}

export function fmt(v, digits = 1) {
  if (!isFinite(v)) return '–';
  return v.toFixed(digits);
}

// ---------------------------------------------------------------------------
// Overview
// ---------------------------------------------------------------------------

export function overviewGeometry(W, H, L, maxDepth) {
  const padL = 46;
  const padR = 14;
  const chartTop = 6;
  const chartBottom = 70;
  const secTop = 82;
  const yWater = secTop + 26;
  const secBottom = H - 26;
  const pxPerM = (secBottom - yWater - 6) / Math.max(1, maxDepth);
  const span = W - padL - padR;
  return {
    padL, padR, chartTop, chartBottom, secTop, yWater, secBottom, pxPerM, span,
    xOf: (xm) => padL + (xm / L) * span,
    xmOf: (px) => ((px - padL) / span) * L,
    yOfDepth: (d) => yWater + d * pxPerM,
  };
}

const LEVEL_EXAG_STEPS = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];

export function drawOverview(ctx, W, H, v, pal) {
  const { L, coarse, depthAt, levelAt, win, dir, U10, ends, editing, points, maxDepth, hoverX } = v;
  const g = overviewGeometry(W, H, L, maxDepth);
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = pal.sky;
  ctx.fillRect(0, 0, W, H);

  // --- top band: Hs(x) and Tp(x) ---------------------------------------
  const M = coarse.M;
  let hsMax = 0;
  let tpMax = 0;
  for (let j = 0; j < M; j++) { hsMax = Math.max(hsMax, coarse.Hs[j]); tpMax = Math.max(tpMax, coarse.Tp[j]); }
  const hsScale = Math.max(0.5, Math.ceil(Math.max(hsMax, v.compare ? v.compare.hsMax : 0) * 2) / 2);
  const tpScale = Math.max(2, Math.ceil(tpMax));
  const bandH = g.chartBottom - g.chartTop;
  const yHs = (h) => g.chartBottom - (h / hsScale) * bandH;
  const yTp = (t) => g.chartBottom - (t / tpScale) * bandH;
  ctx.strokeStyle = pal.grid;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(g.padL, g.chartBottom + 0.5);
  ctx.lineTo(W - g.padR, g.chartBottom + 0.5);
  ctx.stroke();
  // Hs area
  ctx.beginPath();
  ctx.moveTo(g.xOf(0), g.chartBottom);
  for (let j = 0; j < M; j++) ctx.lineTo(g.xOf(j * coarse.dx), yHs(coarse.Hs[j]));
  ctx.lineTo(g.xOf((M - 1) * coarse.dx), g.chartBottom);
  ctx.closePath();
  ctx.fillStyle = pal.accentSoft;
  ctx.fill();
  ctx.beginPath();
  for (let j = 0; j < M; j++) {
    const x = g.xOf(j * coarse.dx);
    const y = yHs(coarse.Hs[j]);
    if (j === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.strokeStyle = pal.accent;
  ctx.lineWidth = 1.5;
  ctx.stroke();
  // breaking marks on the Hs curve
  ctx.fillStyle = pal.breaking;
  for (let j = 0; j < M; j++) if (coarse.breaking[j]) ctx.fillRect(g.xOf(j * coarse.dx) - 0.75, g.chartBottom - 4, 1.5, 4);
  // Tp line
  ctx.beginPath();
  for (let j = 0; j < M; j++) {
    const x = g.xOf(j * coarse.dx);
    const y = yTp(coarse.Tp[j]);
    if (j === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.setLineDash([3, 3]);
  ctx.strokeStyle = pal.muted;
  ctx.lineWidth = 1;
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.font = `11px ${pal.mono}`;
  ctx.textBaseline = 'top';
  ctx.textAlign = 'left';
  ctx.fillStyle = pal.accent;
  ctx.fillText(`Hs  0–${fmt(hsScale, 1)} m   max ${fmt(hsMax, 2)} m`, g.padL + 4, g.chartTop + 1);
  ctx.fillStyle = pal.muted;
  ctx.fillText(`Tp  0–${tpScale} s (dashed)   max ${fmt(tpMax, 1)} s`, g.padL + 4, g.chartTop + 14);
  // comparison lake: Hs along the same fetch distance, dashed
  if (v.compare) {
    const cc = v.compare.coarse;
    const cdir = v.compare.dir;
    ctx.beginPath();
    let started = false;
    for (let j = 0; j < cc.M; j++) {
      const fetch = cdir > 0 ? j * cc.dx : (cc.M - 1 - j) * cc.dx;      // distance from its upwind shore
      const xm = dir > 0 ? fetch : L - fetch;                             // same distance on this lake
      if (xm < 0 || xm > L) continue;
      const x = g.xOf(xm);
      const y = yHs(Math.min(hsScale, cc.Hs[j]));
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    }
    ctx.setLineDash([6, 3]);
    ctx.strokeStyle = pal.breaking;
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = pal.breaking;
    ctx.textAlign = 'left';
    ctx.fillText(`dashed: ${v.compare.name}, same wind, max Hs ${fmt(v.compare.hsMax, 2)} m`, g.padL + 4, g.chartTop + 27);
  }

  // --- section: water, bed, level ---------------------------------------
  const cols = Math.round(g.span);
  const xs = new Float64Array(cols + 1);
  const bedY = new Float64Array(cols + 1);
  const lvlY = new Float64Array(cols + 1);
  let etaMax = 0;
  const etas = new Float64Array(cols + 1);
  for (let c = 0; c <= cols; c++) {
    const xm = (c / cols) * L;
    xs[c] = g.padL + c;
    bedY[c] = g.yOfDepth(depthAt(xm));
    etas[c] = levelAt ? levelAt(xm) : 0;
    etaMax = Math.max(etaMax, Math.abs(etas[c]));
  }
  let E = LEVEL_EXAG_STEPS[0];
  for (const e of LEVEL_EXAG_STEPS) if (etaMax * e * g.pxPerM <= 20) E = e;
  for (let c = 0; c <= cols; c++) lvlY[c] = Math.max(g.chartBottom + 4, Math.min(bedY[c], g.yWater - etas[c] * E * g.pxPerM));

  // water
  const grad = ctx.createLinearGradient(0, g.yWater, 0, g.secBottom);
  grad.addColorStop(0, pal.waterTop);
  grad.addColorStop(1, pal.waterDeep);
  ctx.beginPath();
  ctx.moveTo(xs[0], lvlY[0]);
  for (let c = 1; c <= cols; c++) ctx.lineTo(xs[c], lvlY[c]);
  for (let c = cols; c >= 0; c--) ctx.lineTo(xs[c], bedY[c]);
  ctx.closePath();
  ctx.fillStyle = grad;
  ctx.fill();
  // bed
  ctx.beginPath();
  ctx.moveTo(xs[0], bedY[0]);
  for (let c = 1; c <= cols; c++) ctx.lineTo(xs[c], bedY[c]);
  ctx.lineTo(xs[cols], g.secBottom);
  ctx.lineTo(xs[0], g.secBottom);
  ctx.closePath();
  ctx.fillStyle = pal.bed;
  ctx.fill();
  ctx.beginPath();
  for (let c = 0; c <= cols; c++) { if (c === 0) ctx.moveTo(xs[c], bedY[c]); else ctx.lineTo(xs[c], bedY[c]); }
  ctx.strokeStyle = pal.bedDark;
  ctx.lineWidth = 1.5;
  ctx.stroke();
  // still-water reference and level line
  if (levelAt) {
    ctx.setLineDash([4, 4]);
    ctx.strokeStyle = pal.muted;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(g.padL, g.yWater + 0.5);
    ctx.lineTo(W - g.padR, g.yWater + 0.5);
    ctx.stroke();
    ctx.setLineDash([]);
  }
  ctx.beginPath();
  for (let c = 0; c <= cols; c++) { if (c === 0) ctx.moveTo(xs[c], lvlY[c]); else ctx.lineTo(xs[c], lvlY[c]); }
  ctx.strokeStyle = pal.fg;
  ctx.lineWidth = 1.25;
  ctx.stroke();
  // exposed lakebed: columns where the wind has left only the film on a bed that is normally under water
  {
    let runStart = -1;
    const spans = [];
    for (let c = 0; c <= cols + 1; c++) {
      const xm = (Math.min(c, cols) / cols) * L;
      const still = depthAt(xm);
      const dry = c <= cols && still > 1.5 && still + etas[c] <= 1.0;
      if (dry && runStart < 0) runStart = c;
      if (!dry && runStart >= 0) { spans.push([runStart, c - 1]); runStart = -1; }
    }
    for (const [a, b] of spans) {
      if (b - a < 2) continue;
      ctx.beginPath();
      for (let c = a; c <= b; c++) { if (c === a) ctx.moveTo(xs[c], bedY[c] - 1); else ctx.lineTo(xs[c], bedY[c] - 1); }
      ctx.strokeStyle = pal.foam;
      ctx.globalAlpha = 0.7;
      ctx.lineWidth = 3;
      ctx.stroke();
      ctx.globalAlpha = 1;
      if (b - a > 60) {
        ctx.font = `11px ${pal.mono}`;
        ctx.fillStyle = pal.fg;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText('lakebed exposed by the wind', (xs[a] + xs[b]) / 2, bedY[Math.floor((a + b) / 2)] + 6);
      }
    }
  }

  // breaking zones on the surface
  ctx.strokeStyle = pal.breaking;
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let j = 0; j < M; j++) {
    if (!coarse.breaking[j]) continue;
    const x = g.xOf(j * coarse.dx);
    const c = Math.max(0, Math.min(cols, Math.round(x - g.padL)));
    ctx.moveTo(x, lvlY[c] - 3);
    ctx.lineTo(x, lvlY[c] - 9);
  }
  ctx.stroke();

  // --- axes ---------------------------------------------------------------
  ctx.font = `11px ${pal.mono}`;
  ctx.fillStyle = pal.muted;
  ctx.strokeStyle = pal.grid;
  ctx.lineWidth = 1;
  const dStep = niceStep(maxDepth, 5);
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let d = 0; d <= maxDepth + 1e-9; d += dStep) {
    const y = g.yOfDepth(d);
    ctx.fillText(`${Math.round(d)}`, g.padL - 6, y);
    ctx.beginPath();
    ctx.moveTo(g.padL - 3, y + 0.5);
    ctx.lineTo(g.padL, y + 0.5);
    ctx.stroke();
  }
  ctx.save();
  ctx.translate(12, (g.yWater + g.secBottom) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.fillText('depth, m', 0, 0);
  ctx.restore();
  const kmStep = niceStep(L / 1e3, 8);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (let km = 0; km <= L / 1e3 + 1e-9; km += kmStep) {
    const x = g.xOf(km * 1e3);
    ctx.fillText(`${Math.round(km)}`, x, g.secBottom + 6);
    ctx.beginPath();
    ctx.moveTo(x + 0.5, g.secBottom);
    ctx.lineTo(x + 0.5, g.secBottom + 3);
    ctx.stroke();
  }
  ctx.textAlign = 'right';
  ctx.fillText('km', W - g.padR, g.secBottom + 6);
  // end labels
  ctx.font = `600 12px ${pal.sans}`;
  ctx.fillStyle = pal.fg;
  ctx.textBaseline = 'top';
  ctx.textAlign = 'left';
  ctx.fillText(ends[0], g.padL + 4, g.secTop);
  ctx.textAlign = 'right';
  ctx.fillText(ends[1], W - g.padR - 4, g.secTop);
  // level exaggeration tag
  if (levelAt && E !== 1) {
    ctx.font = `11px ${pal.mono}`;
    ctx.fillStyle = pal.muted;
    ctx.textAlign = 'center';
    ctx.fillText(`water level drawn ×${E}`, g.padL + g.span * (dir > 0 ? 0.72 : 0.28), g.secTop);
  }
  // wind arrow
  if (U10 > 0) {
    const y = g.secTop + 9;
    const len = Math.min(120, 30 + U10 * 3);
    const cx = g.padL + g.span * (dir > 0 ? 0.3 : 0.7);
    const xStart = dir > 0 ? cx - len / 2 : cx + len / 2;
    const xEnd = xStart + dir * len;
    ctx.strokeStyle = pal.wind;
    ctx.fillStyle = pal.wind;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(xStart, y);
    ctx.lineTo(xEnd, y);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(xEnd, y);
    ctx.lineTo(xEnd - dir * 8, y - 4);
    ctx.lineTo(xEnd - dir * 8, y + 4);
    ctx.closePath();
    ctx.fill();
    ctx.font = `11px ${pal.mono}`;
    ctx.textAlign = dir > 0 ? 'right' : 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(`wind ${fmt(U10, 1)} m/s`, xStart - dir * 6, y);
  }

  // --- window marker --------------------------------------------------------
  const wx0 = g.xOf(win.center - win.width / 2);
  const wx1 = g.xOf(win.center + win.width / 2);
  const wW = Math.max(8, wx1 - wx0);
  const wx = (wx0 + wx1) / 2 - wW / 2;
  ctx.fillStyle = pal.accent;
  ctx.globalAlpha = 0.18;
  ctx.fillRect(wx, g.chartTop, wW, g.secBottom - g.chartTop);
  ctx.globalAlpha = 1;
  ctx.strokeStyle = pal.accent;
  ctx.lineWidth = 1.5;
  ctx.strokeRect(wx + 0.5, g.chartTop + 0.5, wW - 1, g.secBottom - g.chartTop - 1);
  ctx.fillRect(wx - 4, g.chartBottom + 2, wW + 8, 8);

  // --- bed control points -------------------------------------------------
  if (editing) {
    ctx.lineWidth = 1.5;
    for (let i = 0; i < points.length; i++) {
      const [xkm, d] = points[i];
      const x = g.xOf(xkm * 1e3);
      const y = g.yOfDepth(d);
      ctx.beginPath();
      ctx.arc(x, y, 5, 0, Math.PI * 2);
      ctx.fillStyle = pal.panel;
      ctx.fill();
      ctx.strokeStyle = pal.bedDark;
      ctx.stroke();
    }
  }

  // --- hover crosshair ------------------------------------------------------
  if (hoverX != null) {
    const x = g.xOf(hoverX);
    ctx.strokeStyle = pal.fg;
    ctx.globalAlpha = 0.5;
    ctx.setLineDash([2, 3]);
    ctx.beginPath();
    ctx.moveTo(x + 0.5, g.chartTop);
    ctx.lineTo(x + 0.5, g.secBottom);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
  }
  return g;
}

// ---------------------------------------------------------------------------
// Zoom window
// ---------------------------------------------------------------------------

/**
 * Window scales.  The vertical scale is the horizontal one times the chosen
 * exaggeration, reduced automatically when the sea is taller than the frame
 * (crests reach about 1.3 Hs with the sharpened crests).
 */
export function zoomGeometry(W, H, width, exag, hsMax = 0, zRef = 0) {
  const sx = W / width;
  const yWater = Math.round(H * 0.4);          // pixel row of the local mean water level
  let sy = sx * exag;
  let auto = false;
  if (hsMax > 0) {
    const fit = (yWater - 46) / (1.3 * hsMax);
    if (fit < sy) { sy = fit; auto = true; }
  }
  // z is measured from the still-water datum; the frame is centred on zRef,
  // the mean level at the window (which wind set-up can move by metres)
  const yOf = (z) => yWater - (z - zRef) * sy;
  return { sx, sy, yWater, exagEff: sy / sx, auto, zRef, yOf };
}

export function fmtExag(g) {
  const e = g.exagEff;
  const txt = e >= 1 ? (Math.abs(e - Math.round(e)) < 0.05 ? String(Math.round(e)) : e.toFixed(1)) : e.toFixed(2);
  return `×${txt}${g.auto ? ' auto-fit' : ''}`;
}

export function drawZoom(ctx, W, H, v, pal) {
  const { fine, eta, win, exag, particles, dir, U10, L, ends, opts, center } = v;
  const n = fine.n;
  const eta0c = fine.eta0[Math.floor(n / 2)];
  const g = zoomGeometry(W, H, win.width, exag, v.hsMax, eta0c);
  const xPx = (p) => (p / (n - 1)) * W;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = pal.sky;
  ctx.fillRect(0, 0, W, H);

  const surfY = new Float64Array(n);
  const bedY = new Float64Array(n);
  let bedMin = Infinity;
  for (let p = 0; p < n; p++) {
    surfY[p] = g.yOf(eta[p]);
    const still = fine.d[p] - fine.eta0[p];
    bedY[p] = g.yOf(-still);
    bedMin = Math.min(bedMin, bedY[p]);
  }

  // water
  const grad = ctx.createLinearGradient(0, g.yWater - 10, 0, H);
  grad.addColorStop(0, pal.waterTop);
  grad.addColorStop(1, pal.waterDeep);
  const waterPath = () => {
    ctx.beginPath();
    ctx.moveTo(0, surfY[0]);
    for (let p = 1; p < n; p++) ctx.lineTo(xPx(p), surfY[p]);
    for (let p = n - 1; p >= 0; p--) ctx.lineTo(xPx(p), Math.min(H + 2, bedY[p]));
    ctx.closePath();
  };
  waterPath();
  ctx.fillStyle = grad;
  ctx.fill();
  // bed
  ctx.beginPath();
  ctx.moveTo(0, Math.min(H + 2, bedY[0]));
  for (let p = 1; p < n; p++) ctx.lineTo(xPx(p), Math.min(H + 2, bedY[p]));
  ctx.lineTo(W, H + 2);
  ctx.lineTo(0, H + 2);
  ctx.closePath();
  ctx.fillStyle = pal.bed;
  ctx.fill();
  ctx.beginPath();
  for (let p = 0; p < n; p++) { const y = Math.min(H + 2, bedY[p]); if (p === 0) ctx.moveTo(0, y); else ctx.lineTo(xPx(p), y); }
  ctx.strokeStyle = pal.bedDark;
  ctx.lineWidth = 2;
  ctx.stroke();

  // HUD text is composed here so that other labels can keep clear of its box
  const pc = Math.floor(n / 2);
  const from = dir > 0 ? ends[0] : ends[1];
  const dist = dir > 0 ? center : L - center;
  const dc = fine.d[pc];
  const stillC = dc - fine.eta0[pc];
  const blownDry = opts.setup && dc <= 1.0 && stillC > 1.5;
  const hud = [
    { text: `${fmt(center / 1e3, 1)} km from ${ends[0]}  ·  fetch ${fmt(dist / 1e3, 0)} km from ${from}`, font: `600 12px ${pal.sans}`, color: pal.fg },
    { text: `depth ${fmt(dc, 1)} m   Hs ${fmt(fine.Hs[pc], 2)} m   Tp ${fmt(fine.Tp[pc], 1)} s   L ${fmt(v.Lp, 0)} m`, font: `11px ${pal.mono}`, color: pal.muted },
  ];
  if (v.Lp > 0 && win.width < 2 * v.Lp && !blownDry) {
    hud.push({ text: `window is ${fmt(win.width / v.Lp, 1)} wavelengths wide: the surface heaves and tilts as each crest passes · F fits 4 waves`, font: `11px ${pal.mono}`, color: pal.wind });
  }
  if (v.vessel) {
    const { def, pose } = v.vessel;
    let line = `${def.name} · draft ${def.draft} m`;
    let color = pal.muted;
    if (pose.aground) {
      line += ` · AGROUND, ${fmt(def.draft + pose.waterUnderKeel, 1)} m of water`;
      color = pal.breaking;
    } else if (v.Lp > 0) {
      line += ` · ${fmt(def.LOA / v.Lp, 1)} wavelengths along the hull · pitch ${fmt((pose.pitch * 180) / Math.PI, 1)}° · heave ${pose.heave >= 0 ? '+' : ''}${fmt(pose.heave, 1)} m`;
    }
    hud.push({ text: line, font: `11px ${pal.mono}`, color });
  }
  // wrap lines that do not fit the canvas at their " · " separators
  const hudLines = [];
  for (const h of hud) {
    ctx.font = h.font;
    if (ctx.measureText(h.text).width <= W - 28) { hudLines.push(h); continue; }
    const parts = h.text.split(' · ');
    let cur = '';
    for (const part of parts) {
      const trial = cur ? `${cur} · ${part}` : part;
      if (cur && ctx.measureText(trial).width > W - 28) { hudLines.push({ ...h, text: cur }); cur = part; } else cur = trial;
    }
    if (cur) hudLines.push({ ...h, text: cur });
  }
  let hudW = 0;
  for (const h of hudLines) { ctx.font = h.font; hudW = Math.max(hudW, ctx.measureText(h.text).width); }
  const hudH = 8 + hudLines.length * 15 + 2;
  const hudBottom = 4 + hudH;

  // local mean level (the frame's reference) and the still-water datum
  ctx.setLineDash([4, 4]);
  ctx.lineWidth = 1;
  ctx.strokeStyle = pal.muted;
  ctx.beginPath();
  ctx.moveTo(0, g.yWater + 0.5);
  ctx.lineTo(W, g.yWater + 0.5);
  ctx.stroke();
  const yStill = g.yOf(0);
  if (opts.setup && Math.abs(yStill - g.yWater) > 2) {
    ctx.strokeStyle = pal.wind;
    ctx.beginPath();
    ctx.moveTo(0, yStill + 0.5);
    ctx.lineTo(W, yStill + 0.5);
    ctx.stroke();
    ctx.setLineDash([]);
    if (yStill > -4 && yStill < H + 4) {
      ctx.font = `11px ${pal.mono}`;
      ctx.fillStyle = pal.wind;
      ctx.textAlign = 'right';
      const above = eta0c < 0;      // datum lies above the water when the wind has lowered it
      const label = `still-water datum, ${fmt(Math.abs(eta0c), 1)} m ${above ? 'above' : 'below'} the water here`;
      let ly = yStill + (above ? -3 : 3);
      ctx.textBaseline = above ? 'bottom' : 'top';
      if (ly < hudBottom + 14 && W - 14 - ctx.measureText(label).width < hudW + 20) { ly = hudBottom + 4; ctx.textBaseline = 'top'; }
      ctx.fillText(label, W - 14, Math.max(2, Math.min(H - 2, ly)));
    }
  }
  ctx.setLineDash([]);

  // particle orbits
  if (opts.orbits && particles) {
    for (const pt of particles) {
      if (pt.trail.length > 1) {
        ctx.beginPath();
        for (let i = 0; i < pt.trail.length; i++) {
          const [tx, ty] = pt.trail[i];
          const px = tx * g.sx;
          const py = g.yOf(ty);
          if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
        }
        ctx.strokeStyle = pal.foam;
        ctx.globalAlpha = 0.65;
        ctx.lineWidth = 1;
        ctx.stroke();
        ctx.globalAlpha = 1;
      }
      ctx.beginPath();
      ctx.arc(pt.cur[0] * g.sx, g.yOf(pt.cur[1]), 2.5, 0, Math.PI * 2);
      ctx.fillStyle = pal.foam;
      ctx.fill();
    }
  }

  // surface line
  ctx.beginPath();
  for (let p = 0; p < n; p++) { if (p === 0) ctx.moveTo(0, surfY[p]); else ctx.lineTo(xPx(p), surfY[p]); }
  ctx.strokeStyle = pal.fg;
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // foam on breaking crests
  ctx.strokeStyle = pal.foam;
  ctx.lineWidth = 3;
  ctx.lineCap = 'round';
  for (let p = 2; p < n - 2; p++) {
    if (!fine.breaking[p]) continue;
    if (eta[p] >= eta[p - 1] && eta[p] >= eta[p + 1] && eta[p] > eta[p - 2] && eta[p] > eta[p + 2]) {
      const span = Math.max(6, Math.min(40, 0.12 * fine.Tp[p] * fine.Tp[p] * g.sx));
      ctx.beginPath();
      ctx.moveTo(xPx(p) - span * (dir > 0 ? 0.2 : 0.8), surfY[p] - 1);
      ctx.lineTo(xPx(p) + span * (dir > 0 ? 0.8 : 0.2), surfY[p] - 1);
      ctx.stroke();
    }
  }
  ctx.lineCap = 'butt';

  // vessel
  if (v.vessel) drawVessel(ctx, W, H, g, v, pal, grad, waterPath);

  // depth ticks (left) and scale bar (bottom right)
  ctx.font = `11px ${pal.mono}`;
  ctx.fillStyle = pal.fg;
  ctx.strokeStyle = pal.grid;
  ctx.lineWidth = 1;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'middle';
  const dStep = niceStep((H - g.yWater) / g.sy, 5);
  for (let d = 0; g.yWater + d * g.sy < H - 58; d += dStep) {
    const y = g.yWater + d * g.sy;
    ctx.beginPath();
    ctx.moveTo(0, y + 0.5);
    ctx.lineTo(8, y + 0.5);
    ctx.stroke();
    ctx.fillText(d === 0 ? '0 m, mean surface' : `-${+d.toFixed(2)} m`, 11, y);
  }
  // one peak wavelength, bracketed above the scale bar
  if (v.Lp > 0 && !blownDry) {
    const y = H - 34;
    const lpPx = v.Lp * g.sx;
    const fits = lpPx <= W - 28;
    const x0 = 14;
    const x1 = fits ? x0 + lpPx : W - 14;
    ctx.strokeStyle = pal.fg;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(x0, y);
    ctx.lineTo(x1, y);
    ctx.moveTo(x0 + 0.5, y - 5);
    ctx.lineTo(x0 + 0.5, y + 5);
    if (fits) { ctx.moveTo(x1 - 0.5, y - 5); ctx.lineTo(x1 - 0.5, y + 5); } else { ctx.moveTo(x1, y); ctx.lineTo(x1 - 7, y - 4); ctx.moveTo(x1, y); ctx.lineTo(x1 - 7, y + 4); }
    ctx.stroke();
    ctx.font = `11px ${pal.mono}`;
    ctx.fillStyle = pal.fg;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    ctx.fillText(fits ? `one peak wavelength, ${fmt(v.Lp, 0)} m` : `one peak wavelength is ${fmt(v.Lp, 0)} m, wider than this window`, x0, y - 6);
  }

  const barM = niceStep(win.width / 4, 1);
  const barPx = barM * g.sx;
  ctx.strokeStyle = pal.fg;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(W - 14 - barPx, H - 14);
  ctx.lineTo(W - 14, H - 14);
  ctx.stroke();
  ctx.textAlign = 'right';
  ctx.textBaseline = 'bottom';
  ctx.fillStyle = pal.fg;
  ctx.fillText(`${+barM.toFixed(2)} m  ·  vertical ${fmtExag(g)}`, W - 14, H - 18);
  if (bedMin > H) {
    const dc = fine.d[Math.floor(n / 2)];
    ctx.textAlign = 'center';
    ctx.fillStyle = pal.foam;
    ctx.fillText(`bed ${fmt(dc, 0)} m down, below this frame ↓`, W / 2, H - 18);
  }

  // blown-dry banner: the bed here is normally under water and the wind has emptied it
  if (blownDry) {
    const lines = [
      'Blown dry',
      `The wind has pushed the water out of this part of the lake: ${fmt(stillC, 1)} m deep in still water, ${fmt(dc, 1)} m film left now.`,
      'Move the window downwind to where the water has piled up, or ease the wind and watch it flow back.',
    ];
    ctx.font = `600 13px ${pal.sans}`;
    let bw = ctx.measureText(lines[0]).width;
    ctx.font = `12px ${pal.sans}`;
    for (const l of lines.slice(1)) bw = Math.max(bw, ctx.measureText(l).width);
    bw = Math.min(W - 24, bw + 24);
    const bh = 70;
    const bx = (W - bw) / 2;
    const by = Math.min(H - bh - 30, Math.max(hudBottom + 10, g.yWater + 24));
    ctx.fillStyle = pal.panel;
    ctx.globalAlpha = 0.92;
    ctx.fillRect(bx, by, bw, bh);
    ctx.globalAlpha = 1;
    ctx.strokeStyle = pal.breaking;
    ctx.lineWidth = 1;
    ctx.strokeRect(bx + 0.5, by + 0.5, bw - 1, bh - 1);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillStyle = pal.breaking;
    ctx.font = `600 13px ${pal.sans}`;
    ctx.fillText(lines[0], bx + 12, by + 10);
    ctx.fillStyle = pal.fg;
    ctx.font = `12px ${pal.sans}`;
    ctx.fillText(lines[1], bx + 12, by + 30, bw - 24);
    ctx.fillStyle = pal.muted;
    ctx.fillText(lines[2], bx + 12, by + 48, bw - 24);
  }

  // HUD, on a translucent backing so it stays readable over a tall ship
  ctx.fillStyle = pal.sky;
  ctx.globalAlpha = 0.82;
  ctx.fillRect(4, 4, hudW + 12, hudH);
  ctx.globalAlpha = 1;
  ctx.textBaseline = 'top';
  ctx.textAlign = 'left';
  hudLines.forEach((h, i) => { ctx.font = h.font; ctx.fillStyle = h.color; ctx.fillText(h.text, 10, 8 + i * 15 + (i > 0 ? 3 : 0)); });
  // wind arrow
  if (U10 > 0) {
    const y = 16;
    const len = Math.min(90, 24 + U10 * 2.2);
    const x0 = dir > 0 ? W - 24 - len : W - 24;
    const x1 = x0 + dir * len;
    ctx.strokeStyle = pal.wind;
    ctx.fillStyle = pal.wind;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(x0, y);
    ctx.lineTo(x1, y);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x1, y);
    ctx.lineTo(x1 - dir * 8, y - 4);
    ctx.lineTo(x1 - dir * 8, y + 4);
    ctx.closePath();
    ctx.fill();
    ctx.textAlign = 'right';
    ctx.fillText(`wind ${fmt(U10, 1)} m/s`, W - 24, y + 8);
  }
  return g;
}

function drawVessel(ctx, W, H, g, v, pal, waterGrad, waterPath) {
  const { def, pose, shapes } = v.vessel;
  const xc = (v.win.center - v.fine.x0);        // metres from the window's left edge
  const cs = Math.cos(pose.pitch);
  const sn = Math.sin(pose.pitch);
  const hs = pose.headingSign;
  // ship coords (x bow+, z up) → window metres → pixels
  const toPx = ([x, z]) => {
    const X = hs * x;
    const Xr = X * cs - z * sn;
    const Zr = X * sn + z * cs;
    return [(xc + Xr) * g.sx, g.yOf(pose.heave + Zr)];
  };
  const path = (pts) => {
    ctx.beginPath();
    pts.forEach((p, i) => { const [px, py] = toPx(p); if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py); });
    ctx.closePath();
  };
  const hull = shapes.find((s) => s.part === 'hull');
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const s of shapes) for (const p of s.pts) { const [px, py] = toPx(p); minX = Math.min(minX, px); maxX = Math.max(maxX, px); minY = Math.min(minY, py); maxY = Math.max(maxY, py); }
  if (maxX < -10 || minX > W + 10) return;
  const fill = { hull: pal.shipHull, boot: pal.shipBoot, house: pal.shipHouse, funnel: pal.shipFunnel, fitting: pal.shipFitting, mast: pal.shipFitting };
  for (const s of shapes) {
    if (s.part === 'deckline') {
      ctx.beginPath();
      const [ax, ay] = toPx(s.pts[0]);
      const [bx, by] = toPx(s.pts[1]);
      ctx.moveTo(ax, ay);
      ctx.lineTo(bx, by);
      ctx.strokeStyle = pal.shipFitting;
      ctx.lineWidth = 1;
      ctx.stroke();
      continue;
    }
    if (s.clip) {
      ctx.save();
      path(hull.pts);
      ctx.clip();
      path(s.pts);
      ctx.fillStyle = fill[s.part];
      ctx.fill();
      ctx.restore();
      continue;
    }
    path(s.pts);
    ctx.fillStyle = fill[s.part] || pal.shipFitting;
    ctx.fill();
    if (s.part === 'hull' || s.part === 'house') {
      ctx.strokeStyle = s.part === 'hull' ? pal.shipEdge : pal.shipFitting;
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }
  // the underwater hull seen through the water: tint the submerged part
  ctx.save();
  ctx.beginPath();
  const x0 = Math.max(0, minX - 2);
  ctx.rect(x0, Math.max(0, minY - 2), Math.min(W, maxX + 2) - x0, Math.min(H, maxY + 2) - Math.max(0, minY - 2));
  ctx.clip();
  waterPath();
  ctx.globalAlpha = 0.6;
  ctx.fillStyle = waterGrad;
  ctx.fill();
  ctx.globalAlpha = 1;
  ctx.restore();
}

// ---------------------------------------------------------------------------
// Spectrum
// ---------------------------------------------------------------------------

export function drawSpectrum(ctx, W, H, v, pal) {
  const { f, sHere, sDeep, fpHere, fpDeep, hsHere, hsDeep } = v;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = pal.sky;
  ctx.fillRect(0, 0, W, H);
  const padL = 34, padR = 28, padT = 10, padB = 26;
  const fMaxPlot = Math.min(f[f.length - 1], Math.max(0.3, 3.2 * Math.max(fpHere || 0, fpDeep || 0)));
  let sMax = 0;
  for (let i = 0; i < f.length; i++) { if (f[i] <= fMaxPlot) sMax = Math.max(sMax, sHere[i], sDeep[i]); }
  if (!(sMax > 0)) {
    ctx.fillStyle = pal.muted;
    ctx.font = `12px ${pal.sans}`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('no wind, no waves', W / 2, H / 2);
    return;
  }
  const xOf = (fr) => padL + (fr / fMaxPlot) * (W - padL - padR);
  const yOf = (s) => H - padB - (s / (sMax * 1.08)) * (H - padT - padB);
  ctx.strokeStyle = pal.grid;
  ctx.lineWidth = 1;
  ctx.fillStyle = pal.muted;
  ctx.font = `10.5px ${pal.mono}`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const fStep = niceStep(fMaxPlot, 5);
  for (let fr = 0; fr <= fMaxPlot + 1e-9; fr += fStep) {
    const x = xOf(fr);
    ctx.beginPath();
    ctx.moveTo(x + 0.5, padT);
    ctx.lineTo(x + 0.5, H - padB);
    ctx.stroke();
    ctx.fillText(fr === 0 ? '0' : `${+fr.toFixed(2)}`, x, H - padB + 4);
    if (fr > 0) ctx.fillText(`${fmt(1 / fr, 1)} s`, x, H - padB + 14);
  }
  ctx.textAlign = 'right';
  ctx.fillText('f, Hz', W - 2, H - padB + 4);
  ctx.fillText('T, s', W - 2, H - padB + 14);
  ctx.save();
  ctx.translate(10, (padT + H - padB) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('S(f), m²/Hz', 0, 0);
  ctx.restore();
  const plot = (s, color, width, dash) => {
    ctx.beginPath();
    let started = false;
    for (let i = 0; i < f.length; i++) {
      if (f[i] > fMaxPlot) break;
      const x = xOf(f[i]);
      const y = yOf(s[i]);
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    }
    ctx.setLineDash(dash);
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.stroke();
    ctx.setLineDash([]);
  };
  // fill under "here"
  ctx.beginPath();
  ctx.moveTo(xOf(f[0]), yOf(0));
  for (let i = 0; i < f.length && f[i] <= fMaxPlot; i++) ctx.lineTo(xOf(f[i]), yOf(sHere[i]));
  ctx.lineTo(xOf(Math.min(fMaxPlot, f[f.length - 1])), yOf(0));
  ctx.closePath();
  ctx.fillStyle = pal.accentSoft;
  ctx.fill();
  plot(sDeep, pal.muted, 1.25, [4, 3]);
  plot(sHere, pal.accent, 2, []);
  ctx.font = `11px ${pal.mono}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillStyle = pal.accent;
  ctx.fillText(`here: Hs ${fmt(hsHere, 2)} m, Tp ${fmt(1 / fpHere, 1)} s`, padL + 6, padT);
  ctx.fillStyle = pal.muted;
  ctx.fillText(`deep: Hs ${fmt(hsDeep, 2)} m, Tp ${fmt(1 / fpDeep, 1)} s`, padL + 6, padT + 13);
}

// ---------------------------------------------------------------------------
// Gauges
// ---------------------------------------------------------------------------

export function drawGauges(ctx, W, H, v, pal) {
  const { samples, now, hours, ends, seicheT } = v;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = pal.sky;
  ctx.fillRect(0, 0, W, H);
  const padL = 40, padR = 54, padT = 10, padB = 22;
  const t0 = now - hours * 3600;
  let m = 0.05;
  for (const s of samples) if (s.t >= t0) m = Math.max(m, Math.abs(s.west), Math.abs(s.east));
  const yMax = m * 1.15;
  const xOf = (t) => padL + ((t - t0) / (hours * 3600)) * (W - padL - padR);
  const yOf = (e) => (padT + H - padB) / 2 - (e / yMax) * ((H - padT - padB) / 2);
  ctx.strokeStyle = pal.grid;
  ctx.lineWidth = 1;
  ctx.fillStyle = pal.muted;
  ctx.font = `10.5px ${pal.mono}`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (let h = -hours; h <= 0; h += hours / 4) {
    const x = xOf(now + h * 3600);
    ctx.beginPath();
    ctx.moveTo(x + 0.5, padT);
    ctx.lineTo(x + 0.5, H - padB);
    ctx.stroke();
    ctx.fillText(h === 0 ? 'now' : `${h} h`, x, H - padB + 5);
  }
  const eStep = niceStep(yMax, 2);
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let e = -Math.floor(yMax / eStep) * eStep; e <= yMax; e += eStep) {
    const y = yOf(e);
    ctx.beginPath();
    ctx.moveTo(padL, y + 0.5);
    ctx.lineTo(W - padR, y + 0.5);
    ctx.strokeStyle = Math.abs(e) < 1e-9 ? pal.muted : pal.grid;
    ctx.stroke();
    ctx.fillText(`${e > 0 ? '+' : ''}${+e.toFixed(2)}`, padL - 5, y);
  }
  ctx.save();
  ctx.translate(10, (padT + H - padB) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.fillText('level, m', 0, 0);
  ctx.restore();
  const line = (key, color) => {
    ctx.beginPath();
    let started = false;
    for (const s of samples) {
      if (s.t < t0) continue;
      const x = xOf(s.t);
      const y = yOf(s[key]);
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    }
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.75;
    ctx.stroke();
    const last = samples[samples.length - 1];
    if (last) {
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(xOf(Math.min(now, last.t)), yOf(last[key]), 3, 0, Math.PI * 2);
      ctx.fill();
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      ctx.font = `11px ${pal.mono}`;
      ctx.fillText(`${last[key] >= 0 ? '+' : ''}${fmt(last[key], 2)}`, W - padR + 6, yOf(last[key]));
    }
  };
  line('west', pal.west);
  line('east', pal.east);
  ctx.fillStyle = pal.muted;
  ctx.font = `11px ${pal.mono}`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`seiche T₁ ≈ ${fmt(seicheT / 3600, 1)} h`, padL + 6, padT);
}
