// bathymetry.js — lake-bed profiles along a cross-section.
//
// Each preset is a long-axis cross-section: distance along the section in km
// and still-water depth in m at that point.  Profiles are drawn through the
// deep parts of each basin, so their mean depth is larger than the lake's
// area-averaged mean depth quoted in the notes.  Values are rounded from NOAA
// Great Lakes bathymetry charts and are meant to be representative, not
// survey-grade.

export const PRESETS = [
  {
    id: 'erie',
    name: 'Lake Erie, Toledo → Buffalo',
    length: 388e3,
    ends: ['Toledo', 'Buffalo'],
    stats: 'mean 19 m, max 64 m',
    note: 'Three basins: a western basin of 7–10 m, a central basin of about 20–25 m, and a small eastern basin that drops to 64 m near Long Point. The prevailing south-westerly wind blows along this axis toward Buffalo, which is why Erie is famous for wind set-up and seiches.',
    points: [[0, 0], [3, 4], [15, 7], [40, 9], [65, 9], [80, 11], [95, 15], [120, 20], [160, 23], [200, 24], [240, 24], [270, 22], [285, 21], [300, 30], [320, 45], [340, 60], [350, 64], [360, 52], [372, 30], [382, 14], [388, 0]],
    // plan-view width (km) along the section: wide western and central basins, narrowing to Buffalo
    widths: [[0, 25], [30, 60], [80, 70], [150, 90], [220, 90], [280, 75], [330, 50], [365, 30], [388, 10]],
  },
  {
    id: 'superior',
    name: 'Lake Superior, Duluth → Whitefish Bay',
    length: 563e3,
    ends: ['Duluth', 'Whitefish Bay'],
    stats: 'mean 147 m, max 406 m',
    note: 'A deep rock basin. Waves generated over its 560 km fetch never feel the bottom until the last few kilometres, and the great depth makes wind set-up tiny for the same wind.',
    points: [[0, 0], [10, 60], [40, 150], [90, 220], [150, 280], [200, 330], [260, 400], [300, 380], [350, 300], [400, 250], [450, 220], [500, 160], [540, 80], [560, 20], [563, 0]],
    widths: [[0, 20], [60, 120], [200, 200], [350, 250], [480, 150], [563, 40]],
  },
  {
    id: 'ontario',
    name: 'Lake Ontario, Hamilton → Kingston',
    length: 311e3,
    ends: ['Hamilton', 'Kingston'],
    stats: 'mean 86 m, max 244 m',
    note: 'Smaller in area than Erie but far deeper: a single bowl that reaches 244 m.',
    points: [[0, 0], [8, 40], [30, 90], [60, 130], [100, 170], [140, 200], [180, 230], [210, 244], [240, 200], [270, 120], [295, 60], [311, 0]],
    widths: [[0, 20], [60, 70], [150, 85], [250, 70], [311, 25]],
  },
  {
    id: 'michigan',
    name: 'Lake Michigan, Mackinac → Chicago',
    length: 494e3,
    ends: ['Mackinac', 'Chicago'],
    stats: 'mean 85 m, max 281 m',
    note: 'A long north–south basin: deep Chippewa basin in the north, a shallower southern basin, and a long gentle shelf toward Chicago.',
    points: [[0, 0], [20, 60], [60, 150], [100, 220], [140, 281], [180, 250], [220, 180], [260, 120], [300, 140], [340, 120], [380, 100], [430, 80], [470, 40], [494, 0]],
    widths: [[0, 60], [100, 130], [250, 120], [400, 140], [494, 40]],
  },
  {
    id: 'flat-shallow',
    name: 'Flat basin, 10 m deep',
    length: 200e3,
    ends: ['West', 'East'],
    stats: 'uniform 10 m',
    note: 'A control case: every wave of more than about 2 s period feels this bottom, so growth, speed and set-up are all depth-limited.',
    points: [[0, 0], [1.5, 10], [198.5, 10], [200, 0]],
  },
  {
    id: 'flat-deep',
    name: 'Flat basin, 200 m deep',
    length: 200e3,
    ends: ['West', 'East'],
    stats: 'uniform 200 m',
    note: 'The deep-water control case with the same length and fetch: compare wave height, period and set-up with the 10 m basin under the same wind.',
    points: [[0, 0], [1.5, 200], [198.5, 200], [200, 0]],
  },
  {
    id: 'beach',
    name: 'Deep basin with a gentle beach',
    length: 100e3,
    ends: ['West', 'East'],
    stats: '60 m, 1:170 slope on the east shore',
    note: 'Waves generated over 90 km of 60 m water run up a long, uniform slope on the downwind side: the classic shoaling-and-breaking picture.',
    points: [[0, 0], [1, 60], [90, 60], [100, 0]],
  },
];

export function presetById(id) {
  return PRESETS.find((p) => p.id === id) || PRESETS[0];
}

/**
 * Shape-preserving piecewise-cubic Hermite interpolant (Fritsch–Carlson).
 * Monotone between control points, so the bed never overshoots the points
 * the way a plain cubic spline would.  xs must be strictly increasing.
 */
export function pchip(xs, ys) {
  const n = xs.length;
  if (n === 1) return () => ys[0];
  const h = new Float64Array(n - 1);
  const delta = new Float64Array(n - 1);
  for (let i = 0; i < n - 1; i++) {
    h[i] = xs[i + 1] - xs[i];
    delta[i] = (ys[i + 1] - ys[i]) / h[i];
  }
  const m = new Float64Array(n);
  if (n === 2) {
    m[0] = m[1] = delta[0];
  } else {
    for (let i = 1; i < n - 1; i++) {
      if (delta[i - 1] * delta[i] <= 0) {
        m[i] = 0;
      } else {
        const w1 = 2 * h[i] + h[i - 1];
        const w2 = h[i] + 2 * h[i - 1];
        m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
      }
    }
    m[0] = endSlope(h[0], h[1], delta[0], delta[1]);
    m[n - 1] = endSlope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
  }
  return function evaluate(x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    let lo = 0;
    let hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (xs[mid] <= x) lo = mid; else hi = mid;
    }
    const t = (x - xs[lo]) / h[lo];
    const t2 = t * t;
    const t3 = t2 * t;
    const h00 = 2 * t3 - 3 * t2 + 1;
    const h10 = t3 - 2 * t2 + t;
    const h01 = -2 * t3 + 3 * t2;
    const h11 = t3 - t2;
    return h00 * ys[lo] + h10 * h[lo] * m[lo] + h01 * ys[lo + 1] + h11 * h[lo] * m[lo + 1];
  };
}

function endSlope(h0, h1, d0, d1) {
  let m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
  if (Math.sign(m) !== Math.sign(d0)) m = 0;
  else if (Math.sign(d0) !== Math.sign(d1) && Math.abs(m) > Math.abs(3 * d0)) m = 3 * d0;
  return m;
}

/** A depth function d(x) in metres for x in metres along the section. */
export function depthFunction(profile) {
  const xs = profile.points.map((p) => p[0] * 1e3);
  const ys = profile.points.map((p) => p[1]);
  const f = pchip(xs, ys);
  return (x) => Math.max(0, f(x));
}

/** Plan-view width b(x) in metres (uniform 1 m channel when the profile has none). */
export function widthFunction(profile) {
  if (!profile.widths) return () => 1;
  const xs = profile.widths.map((p) => p[0] * 1e3);
  const ys = profile.widths.map((p) => p[1] * 1e3);
  const f = pchip(xs, ys);
  return (x) => Math.max(1, f(x));
}

/** n depths sampled at x_j = j/(n−1)·L (cell *nodes*). */
export function sampleDepths(profile, n) {
  const f = depthFunction(profile);
  const out = new Float64Array(n);
  const L = profile.length;
  for (let j = 0; j < n; j++) out[j] = f((j / (n - 1)) * L);
  return out;
}

/** Deep copy of a preset for editing. */
export function cloneProfile(p) {
  return { ...p, points: p.points.map((q) => [q[0], q[1]]), widths: p.widths ? p.widths.map((q) => [q[0], q[1]]) : undefined };
}

/** Rescale a profile's horizontal extent to a new length (m), keeping its shape. */
export function withLength(profile, length) {
  const s = length / profile.length;
  return {
    ...profile,
    length,
    points: profile.points.map((q) => [q[0] * s, q[1]]),
    widths: profile.widths ? profile.widths.map((q) => [q[0] * s, q[1]]) : undefined,
  };
}
