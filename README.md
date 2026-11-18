# Lake Wave Lab

An interactive, physics-based wave simulator for a lake cross-section. Pick a
lake-bed profile (Lake Erie, Superior, Ontario, Michigan, flat control basins,
a beach, or draw your own), set a wind, and watch:

* the sea grow along the fetch and stop growing where the bed limits it,
* each wave component shorten, slow down and shoal as it feels the bottom,
* breaking where the water gets too shallow,
* water particles tracing circles in deep water and flattened ellipses in shallow water,
* the whole lake tilt under the wind (set-up) and slosh back as a seiche when the wind changes.

No build step, no dependencies: plain HTML, CSS and JavaScript on a 2-D canvas.

## Run it

* **Standalone file:** open `dist/wave-sim.html` in any browser (works from disk).
* **Development:** serve the repository root and open `index.html`, which loads the
  ES modules in `src/`:

  ```bash
  npm run serve
  ```

  then visit <http://localhost:8080/>.

* **Rebuild the standalone file** after editing `src/` or `index.html`:

  ```bash
  npm run build
  ```

* **Tests** (Node 20+, no packages to install):

  ```bash
  npm test
  ```

## How it is put together

| File | Role |
| --- | --- |
| `src/physics.js` | Pure formulas with their sources: dispersion relation, group speed, shoaling, Stokes crest term, orbital motion, wind stress, SMB growth curves, duration limit, JONSWAP/TMA spectrum, breaking limits, equilibrium set-up, seiche period. |
| `src/bathymetry.js` | Lake profiles (depth and plan-view width along the long axis) and a shape-preserving PCHIP interpolator. |
| `src/wavefield.js` | Marches a 64-component spectrum across the basin: fetch/depth-limited generation, shoaling, breaking with permanent energy loss, phase integration, and the fine-resolution surface for the zoom window. |
| `src/slowmodel.js` | One-dimensional shallow-water model (staggered grid, forward–backward scheme, quadratic bottom friction, variable width) for wind set-up and seiches. |
| `src/render.js` | Canvas drawing: whole-lake overview, to-scale window, spectrum and gauge charts. |
| `src/main.js` | State, the two clocks, controls, readouts and the generated explanation. |
| `build.mjs` | Inlines the modules into `dist/wave-sim.html`. |
| `test/` | Unit tests that pin the formulas to hand calculations and textbook limits. |

### Two clocks

Wind waves have periods of seconds; set-up and seiches take hours. The window
animates the waves in real time while the lake clock runs 60 to 3600 times
faster for the basin-scale response and for wave growth with wind duration.
Pause the lake clock and the sea becomes quasi-steady (instant response to the
wind slider).

### The physics, layer by layer

1. **Wind to waves.** Wind stress factor U_A = 0.71 U10^1.23 drives the
   Sverdrup–Munk–Bretschneider shallow-water growth curves (Shore Protection
   Manual 1984, eqs. 3-39/3-40) with the mean depth between the upwind shore
   and the point as the generation depth. Duration limit from SPM eq. 3-35.
2. **Spectrum.** JONSWAP (γ = 3.3) × TMA depth factor, normalised to the
   SMB significant height, over 64 log-spaced components with fixed random phases.
3. **Propagation.** Linear (Airy) theory per component: ω² = g k tanh(kd),
   C_g = n C, shoaling K_s = √(C_g,gen / C_g,local).
4. **Breaking.** H_s capped at 0.6 d (individual waves break near 0.78 d) and
   at the Miche steepness limit; the loss is permanent, with a 20 km regrowth
   fetch so the wind can rebuild the sea past a sill.
5. **Crest shape.** Second-order Stokes term per component, capped where the
   Ursell number invalidates Stokes theory.
6. **Orbits.** Linear-theory particle displacements summed over components.
7. **Set-up and seiche.** ∂η/∂t = −(1/b) ∂(b h u)/∂x and
   ∂u/∂t = −g ∂η/∂x + τ/(ρh) − C_f |u| u / h, with τ = ρ_air C_D U10² (Wu 1982).
   The seiche period is a result of the model; the readout also shows the
   long-wave estimate T₁ = 2 ∫ dx / √(g h).

The "Why it looks like this" panel on the page turns the current numbers into
sentences, and the collapsible notes at the bottom list every formula with its source.

### Known limits

* One horizontal dimension: no refraction, directional spreading or islands.
  The funnelling that concentrates storm surge at Buffalo is a 2-D effect.
* The SMB curves are empirical fits from the 1980s; the Coastal Engineering
  Manual prefers friction-velocity scaling for deep water.
* Profiles are rounded from NOAA Great Lakes bathymetry along the long axis,
  so they run deeper than each lake's area-mean depth.

## Towards a 3-D version

The physics modules know nothing about the canvas. A plan-view or 3-D
visualiser would reuse `physics.js` unchanged and replace `wavefield.js` with a
2-D grid version (adding refraction through a ray or wave-action solver and
directional spreading), rendered with WebGL from the same component tables.
