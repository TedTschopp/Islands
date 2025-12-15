# Procedural Island (local copy)

This repo is a self-contained, browser-run procedural landscape + drainage + erosion simulation.
It renders a Voronoi/Delaunay mesh in WebGL, computes a drainage graph each step, then applies simple fluvial erosion/deposition and a small amount of diffusion.

## Run locally

From this folder:

```bash
python3 -m http.server 5173
```

Then open:

- http://localhost:5173/

(You can also open `index.html` directly, but serving via HTTP avoids browser restrictions around some asset loading.)

## File map (where things live)

- `index.html`
  - Main entrypoint.
  - Defines **simulation state arrays** (`bedrock_thickness`, `sediment_thickness`, `lake_thickness`, …).
  - Defines **tunable parameters** (sea level, rainfall, erosion coefficients, timestep, etc.).
  - Orchestrates the per-step workflow in `step(dt)` and schedules animation.
- `voronoi.js`
  - Builds the **Voronoi tiling** on a regular grid of "cells".
  - Produces global arrays like `cell_x`, `cell_y`, `cell_neighbours`, and triangle lists used for rendering.
- `drainage.js`
  - Builds a **directed drainage graph** over cells (including lake filling) and then performs **routing + erosion/deposition** in `run_drainage(dt)`.
  - Produces `discharge[...]` used by river rendering.
- `webgl.js`
  - WebGL setup and all drawing.
  - Draws: terrain triangles, coastlines (lake/sea edges), and rivers (from `discharge`).
- `utilities.js`
  - Small helpers + noise wrapper (`SNoise`) + simple data structures used by drainage.
- `noise.js`
  - External simplex/perlin noise implementation (keep as-is unless you know what you’re doing).

## Coordinate system + units

- Horizontal coordinates:
  - `cell_x[cell]`, `cell_y[cell]` are in **pixel-space** coordinates on the canvas.
  - `pixel_length` (in `index.html`) converts pixel distance to **meters**.
- Vertical coordinates:
  - `bedrock_thickness`, `sediment_thickness`, `lake_thickness`, `cell_elevation` are treated as **meters**.
- Time:
  - `dt` passed to `step(dt)` and `run_drainage(dt)` is **years**.

## The main workflow (where to add new steps)

The per-frame simulation loop is:

1. Update `time` and (optionally) sea level.
2. Apply uplift to bedrock (`uplift(dt)`).
3. Recompute drainage inputs (`initialize_drainage()`).
4. Rebuild the drainage graph (`calculate_drainage_graph()`).
5. Route water/sediment and apply erosion/deposition (`run_drainage(dt)`).
6. Render (`render()`).

This is implemented in `step(dt)` inside `index.html`.

### Adding a new process step

Most new “geologic process” steps belong in one of these places:

- **Before drainage is computed** (changes topography that affects flow):
  - Add just before `initialize_drainage()` in `step(dt)`.
  - Example: tectonic tilting, localized uplift/subsidence, crater impacts.
- **After drainage routing** (post-processing on the updated terrain):
  - Add just after `run_drainage(dt)`.
  - Example: smoothing, coastal processes, post-erosion mass wasting.
- **Inside drainage** (if it depends on discharge/graph structure):
  - Extend `run_drainage(dt)` in `drainage.js`.

Rule of thumb: if your new process changes elevation, you usually want it to happen *before* the next `calculate_drainage_graph()` so the flow field responds.

## “Settings” you’ll most likely tweak

These are declared near the top of `index.html`:

- `pixel_length` (meters): horizontal scale per Voronoi pixel.
- `mean_sea_level` / `sea_level` (meters): vertical datum and current sea level.
- `rainfall` (meters/year): baseline water input per cell area.
- `timestep` (years): how large each `step(dt)` is when you press 1/2/3 or space.
- Erosion/deposition tuning:
  - `fluvial_transport_coefficient`
  - `sedimentation_distance`
  - `sediment_density`
  - `landscape_diffusion_coefficient`
- Display tuning:
  - `lake_threshold` (meters): minimum water depth shown as water.

If you’re changing the look (not the physics):

- Color ramp / shading lives in `webgl.js` fragment shader (`f`).
- River visibility depends on `draw_rivers(lower_threshold, reference, line_width)` in `render()`.

## Debugging tips

- Open DevTools Console.
- If you change any of the global arrays, re-run `initialize_drainage()` + `calculate_drainage_graph()` to rebuild flow.
- If rendering looks wrong, confirm `triangle_count`, `triangle_vertices`, and `triangle_cells` exist (they’re created by `create_voronoi_tiling(...)`).
