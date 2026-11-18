// vessel.js — ships drawn on the surface for scale.
//
// Dimensions are loaded-condition figures typical of each class:
//   tanker:     Great Lakes product tankers (Algoma Hansa / Algoterra size)
//   seawaymax:  the largest hull that fits the St. Lawrence Seaway locks
//   thousand:   the thirteen 1000-footers confined to the upper four lakes
//   classic:    pilothouse-forward straight-decker of the Edmund Fitzgerald era
//   boat:       an 8 m sport-fishing boat, for contrast
//
// The vessel is hydrostatic: it floats at the mean level and mean slope of
// the water along its hull (a least-squares line through the surface under
// it), with no inertia or added mass.  That is exact for a hull much shorter
// or much longer than the waves and a fair sketch in between.

export const VESSELS = [
  { id: 'none', name: 'none' },
  { id: 'tanker', name: 'Great Lakes tanker, 150 m', LOA: 150, beam: 23.7, draft: 8.0, hullDepth: 12.0, layout: 'tanker' },
  { id: 'seawaymax', name: 'Seawaymax laker, 225 m', LOA: 225.5, beam: 23.8, draft: 8.1, hullDepth: 14.3, layout: 'laker-aft' },
  { id: 'thousand', name: '1000-foot laker, 309 m', LOA: 308.7, beam: 32.0, draft: 8.8, hullDepth: 17.0, layout: 'laker-aft' },
  { id: 'classic', name: 'Classic laker, pilothouse forward, 222 m', LOA: 222.5, beam: 22.9, draft: 7.6, hullDepth: 11.9, layout: 'laker-forward' },
  { id: 'boat', name: 'Fishing boat, 8 m', LOA: 8, beam: 2.6, draft: 0.7, hullDepth: 1.6, layout: 'boat' },
];

export function vesselById(id) {
  return VESSELS.find((v) => v.id === id) || VESSELS[0];
}

/**
 * Least-squares line η = b + m·(x − xc) through the surface samples under the
 * hull.  eta/x: fine-grid arrays; pc: centre index; halfLen: half the hull
 * length in metres.  Returns { b, m, samples }.
 */
export function fitWaterline(eta, x, pc, halfLen) {
  const n = eta.length;
  const dx = n > 1 ? x[1] - x[0] : 1;
  const hp = Math.max(0, Math.round(halfLen / dx));
  const p0 = Math.max(0, pc - hp);
  const p1 = Math.min(n - 1, pc + hp);
  let sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  const xc = x[pc];
  for (let p = p0; p <= p1; p++) {
    const u = x[p] - xc;
    sw += 1; sx += u; sy += eta[p]; sxx += u * u; sxy += u * eta[p];
  }
  if (sw < 2) return { b: eta[pc], m: 0, samples: sw };
  const den = sw * sxx - sx * sx;
  const m = den > 1e-12 ? (sw * sxy - sx * sy) / den : 0;
  const b = (sy - m * sx) / sw;
  return { b, m, samples: sw };
}

/**
 * Pose of a vessel centred on fine-grid index pc.
 *   headingSign: +1 bow toward +x, −1 bow toward −x
 * Returns { heave, pitch (rad), aground, waterUnderKeel }.
 */
export function vesselPose(def, fine, eta, pc, headingSign) {
  const fit = fitWaterline(eta, fine.x, pc, def.LOA / 2);
  const stillDepth = fine.d[pc] - fine.eta0[pc];     // bed below still-water level
  const waterUnderKeel = stillDepth + fit.b - def.draft;
  if (waterUnderKeel < 0) {
    // keel on the bed: lift the hull so the keel sits at the bed, no pitch
    return { heave: def.draft - stillDepth, pitch: 0, aground: true, waterUnderKeel, headingSign };
  }
  return { heave: fit.b, pitch: Math.atan(fit.m), aground: false, waterUnderKeel, headingSign };
}

/**
 * Side-view outline in ship coordinates: bow at +x, z = 0 at the loaded
 * waterline, metres.  Returns shapes drawn in order: { pts, part, clip }.
 *   part: 'hull' | 'boot' | 'house' | 'funnel' | 'fitting' | 'mast'
 *   clip: true → clipped to the hull outline (for the boot topping)
 */
export function vesselShapes(def) {
  const L = def.LOA;
  const T = def.draft;
  const F = def.hullDepth - def.draft;          // freeboard
  const shapes = [];
  const rect = (x0, x1, z0, z1, part) => shapes.push({ part, pts: [[x0, z0], [x1, z0], [x1, z1], [x0, z1]] });

  if (def.layout === 'boat') {
    shapes.push({ part: 'hull', pts: [[-L / 2, F], [L / 2, F + 0.35], [L / 2 - 0.6, -0.15], [L / 2 - 2.2, -T], [-L / 2 + 0.6, -T], [-L / 2, -0.15]] });
    shapes.push({ part: 'boot', clip: true, pts: [[-L, -0.25], [L, -0.25], [L, 0.15], [-L, 0.15]] });
    rect(-1.6, 1.4, F, F + 1.5, 'house');
    rect(-0.4, 0.4, F + 1.5, F + 3.2, 'mast');
    return shapes;
  }

  // hull with a raked stem and a cut-away stern
  shapes.push({ part: 'hull', pts: [[-L / 2, F], [L / 2, F], [L / 2 - 0.045 * L, -T], [-L / 2 + 0.07 * L, -T], [-L / 2 + 0.015 * L, -0.35 * T], [-L / 2, 0.2 * F]] });
  shapes.push({ part: 'boot', clip: true, pts: [[-L, -1.2], [L, -1.2], [L, 0.6], [-L, 0.6]] });
  // forecastle
  rect(L / 2 - 0.07 * L, L / 2 - 0.004 * L, F, F + 2.2, 'house');
  rect(L / 2 - 0.035 * L, L / 2 - 0.03 * L, F + 2.2, F + 13, 'mast');

  const aftHouse = (height, decks) => {
    const x0 = -L / 2 + 0.035 * L;
    const x1 = -L / 2 + 0.135 * L;
    rect(x0, x1, F, F + height, 'house');
    rect(x0 + 0.012 * L, x1 - 0.012 * L, F + height, F + height + 3, 'house');
    rect(x0 + 0.005 * L, x0 + 0.03 * L, F + height + 2, F + height + 7, 'funnel');
    rect(x1 - 0.02 * L, x1 - 0.015 * L, F + height + 3, F + height + 10, 'mast');
    for (let i = 1; i < decks; i++) shapes.push({ part: 'deckline', pts: [[x0, F + (height * i) / decks], [x1, F + (height * i) / decks]] });
  };

  if (def.layout === 'tanker') {
    aftHouse(11, 4);
    for (const xf of [-0.18, -0.06, 0.06, 0.18]) rect(xf * L - 0.012 * L, xf * L + 0.012 * L, F, F + 1.6, 'fitting');
    // raised catwalk along the deck
    shapes.push({ part: 'deckline', pts: [[-L / 2 + 0.14 * L, F + 1.8], [L / 2 - 0.08 * L, F + 1.8]] });
    rect(-0.3 * L, -0.29 * L, F, F + 1.8, 'mast');
    rect(0.3 * L, 0.31 * L, F, F + 1.8, 'mast');
  } else if (def.layout === 'laker-aft') {
    aftHouse(12, 4);
    // hatch covers
    const hx0 = -L / 2 + 0.16 * L;
    const hx1 = L / 2 - 0.1 * L;
    const pitch = 0.05 * L;
    for (let x = hx0; x + 0.8 * pitch <= hx1; x += pitch) rect(x, x + 0.8 * pitch, F, F + 0.7, 'fitting');
    // self-unloading boom stowed along the deck, pivoting on a pedestal by the house
    const px = -L / 2 + 0.15 * L;
    rect(px, px + 0.025 * L, F, F + 6, 'fitting');
    shapes.push({ part: 'fitting', pts: [[px + 0.01 * L, F + 6], [px + 0.37 * L, F + 4.6], [px + 0.37 * L, F + 5.6], [px + 0.01 * L, F + 7.4]] });
  } else if (def.layout === 'laker-forward') {
    // pilothouse forward, engine house and funnel aft, hatches between
    const bx0 = L / 2 - 0.12 * L;
    const bx1 = L / 2 - 0.025 * L;
    rect(bx0, bx1, F, F + 6, 'house');
    rect(bx0 + 0.01 * L, bx1 - 0.01 * L, F + 6, F + 9, 'house');
    rect(bx0 + 0.025 * L, bx1 - 0.03 * L, F + 9, F + 11.5, 'house');
    shapes.push({ part: 'deckline', pts: [[bx0, F + 3], [bx1, F + 3]] });
    const ax0 = -L / 2 + 0.03 * L;
    const ax1 = -L / 2 + 0.13 * L;
    rect(ax0, ax1, F, F + 7, 'house');
    rect(ax0 + 0.02 * L, ax0 + 0.045 * L, F + 7, F + 13.5, 'funnel');
    shapes.push({ part: 'deckline', pts: [[ax0, F + 3.5], [ax1, F + 3.5]] });
    rect(ax1 - 0.01 * L, ax1 - 0.005 * L, F + 7, F + 15, 'mast');
    const hx0 = ax1 + 0.03 * L;
    const hx1 = bx0 - 0.02 * L;
    const pitch = 0.045 * L;
    for (let x = hx0; x + 0.75 * pitch <= hx1; x += pitch) rect(x, x + 0.75 * pitch, F, F + 0.6, 'fitting');
  }
  return shapes;
}
