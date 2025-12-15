"use strict";

/*
Drainage graph + fluvial erosion/deposition

This file does two main jobs:

1) Build a directed drainage graph (with lake filling)
    - `initialize_drainage()` resets per-step derived fields.
    - `calculate_drainage_graph()` assigns 1–2 downstream neighbours per cell
       via `best_neighbours[cell*2 + {0,1}]`.

2) Route water/sediment and modify terrain
    - `run_drainage(dt)` processes cells in a topological order derived from
       `in_degree` and applies:
          - rainfall accumulation
          - sediment transport capacity
          - erosion (when capacity > sediment)
          - deposition (when capacity < sediment)
          - a small, targeted "diffusion" term (coastal/underwater smoothing)

Important arrays (all globals shared across files):

- Inputs (from `index.html`):
   - `bedrock_thickness[cell]`, `sediment_thickness[cell]`
   - `sea_level`, `rainfall`, `pixel_length`, tuning coefficients
   - `erodibility[cell]`

- Derived per-step fields:
   - `cell_elevation[cell] = bedrock_thickness + sediment_thickness`
   - `lake_thickness[cell]` (computed by lake filling)
   - `best_neighbours[cell*2 + 0]`: primary downstream neighbour cell id
   - `best_neighbours[cell*2 + 1]`: optional secondary neighbour (braiding/splits)
   - `best_positions[...]`: neighbour-list positions that correspond to those cells
   - `discharge[cell*max_neighbours + i]`: water routed along a specific edge

Where to tweak “settings”:
- Most knobs are defined in `index.html` and used here.
- If you add a new process, `run_drainage(dt)` is the right place when the
   process depends on discharge/gradients/flow direction.
*/

let cell_state = new Uint8Array(cell_count);
const UNTOUCHED = 0;
const LOWER_ADJACENT = 1;
const UPPER_ADJACENT = 2;
const COVERED = 3;
let upper_adjacent_cells = new PriorityQueue(cell_count, cell_elevation);
let lower_adjacent_cells = new RandomQueue(cell_count);
let ready_to_drain_cells = new Stack(cell_count);
let best_neighbours = new Int32Array(cell_count*2);
let best_positions = new Int32Array(cell_count*2);

let in_degree = new Uint8Array(cell_count);
let lake_index = new Int32Array(cell_count);
let next_lake_index = 0;
let collected_water = new Float32Array(cell_count);
let collected_sediment =  new Float32Array(cell_count);
let discharge = new Float32Array(cell_count * max_neighbours);
let covered_time = new Int32Array(cell_count);
let covered_cell_count;

function initialize_drainage()
{
   // Resets all per-step derived fields.
   // Call this after *any* terrain modification if you need fresh drainage.
   upper_adjacent_cells.clear();
   lower_adjacent_cells.clear();
   covered_cell_count = 0;
   next_lake_index = 0;
   for (let i=0; i<cell_count; i++)
   {
      cell_elevation[i] = bedrock_thickness[i]+sediment_thickness[i];
      cell_state[i] = UNTOUCHED;
      lake_thickness[i] = 0;
      in_degree[i] = 0;
      collected_water[i] = 0;
      collected_sediment[i] = 0;
   }
   for (let i=0; i<cell_count*max_neighbours; i++) { discharge[i] = 0; }
   for (let i=0; i<cell_count*2; i++) { best_neighbours[i] = -1; best_positions[i] = 0;}
   for (let i=0; i<cell_width; i++)
   {
      mark_edge_cell(i);
      mark_edge_cell(cell_count - cell_width + i);
   }
   for (let j=0; j<cell_height; j++)
   {
      mark_edge_cell(j*cell_width); // left edge
      mark_edge_cell(j*cell_width + cell_width -1); //right edge
   }
}

function mark_edge_cell(cell)
{
   if (morton_order) { cell = interleave_bits(cell, morton_bits); }
   if (cell_elevation[cell] > sea_level) return;
   if (cell_state[cell] === UNTOUCHED)
   {
      cell_state[cell] = LOWER_ADJACENT;
      lower_adjacent_cells.push(cell);
   }
}


function calculate_drainage_graph()
{
   // Builds the directed drainage relation by "covering" cells from low to high.
   // This behaves like a priority-flood / depression filling pass:
   // - water can exit at edge cells that are below sea level
   // - inland basins become lakes with `lake_thickness`
   let cover_level = sea_level;
   let finished = false;
   while (!finished)
   {
      if (lower_adjacent_cells.is_nonempty())
      {
         let next_cell = lower_adjacent_cells.pop();
         cover_lake(next_cell, cover_level);
         lake_index[next_cell] = next_lake_index;
      }
      else if (upper_adjacent_cells.is_nonempty())
      {
         let next_cell = upper_adjacent_cells.pop();
         next_lake_index = next_cell;
         cover_level = cell_elevation[next_cell];
         cover_non_lake(next_cell, cover_level);
         lake_index[next_cell] = -1;
      }
      else { finished = true; }
   }
}

function cover_non_lake(cell, cover_level)
{
   // Called when we are raising the cover level over dry land.
   // Chooses 1–2 already-covered neighbours as downstream targets.
   cell_state[cell] = COVERED;
   covered_time[cell] = covered_cell_count;
   covered_cell_count += 1;
   let count = cell_neighbour_count[cell];
   let lowest_neighbour = -1;
   let lowest_level = Number.MAX_VALUE;
   let second_lowest_neighbour = -1;
   let second_lowest_level = Number.MAX_VALUE;
   for (let i=0; i<count; i++)
   {
      let position = cell*max_neighbours + i;
      let neighbour = cell_neighbours[position];
      if (cell_state[neighbour] === UNTOUCHED)
      {
         if (cell_elevation[neighbour] <= cover_level)
         {
            cell_state[neighbour] = LOWER_ADJACENT;
            lower_adjacent_cells.push(neighbour);
         }
         else
         {
            cell_state[neighbour] = UPPER_ADJACENT;
            upper_adjacent_cells.push(neighbour);
         }
      }
      else if (cell_state[neighbour] === COVERED)
      {
         let neighbour_level = cell_elevation[neighbour] + lake_thickness[neighbour] - 2*Math.random();
         if (neighbour_level < lowest_level)
         {
            second_lowest_level = lowest_level;
            second_lowest_neighbour = lowest_neighbour;
            best_positions[cell*2+1] = best_positions[cell*2];
            lowest_level = neighbour_level;
            lowest_neighbour = neighbour;
            best_positions[cell*2] = position;
         }
         else if (neighbour_level < second_lowest_level)
         {
            second_lowest_level = neighbour_level;
            second_lowest_neighbour = neighbour;
            best_positions[cell*2+1] = position;
         }
      }
   }
   if (lowest_neighbour != -1)
   {
      best_neighbours[cell*2] = lowest_neighbour;
      in_degree[lowest_neighbour]+=1;
   }
   if (second_lowest_neighbour != -1)
   {
      best_neighbours[cell*2+1] = second_lowest_neighbour;
      in_degree[second_lowest_neighbour]+=1;
   }
}

function cover_lake(cell, cover_level)
{
   // Called when `cell_elevation[cell] <= cover_level` in a filled basin.
   // Sets `lake_thickness[cell]` so the lake surface sits at `cover_level`.
   cell_state[cell] = COVERED;
   covered_time[cell] = covered_cell_count;
   covered_cell_count += 1;
   lake_thickness[cell] = cover_level - cell_elevation[cell];
   let count = cell_neighbour_count[cell];
   let lowest_neighbour = -1;
   let lowest_elevation = Number.MAX_VALUE;
   let last_covered_neighbour = -1;
   let best_covered_time = -1;
   for (let i=0; i<count; i++)
   {
      let position = cell*max_neighbours + i;
      let neighbour = cell_neighbours[position];
      if (cell_state[neighbour] === UNTOUCHED)
      {
         if (cell_elevation[neighbour] <= cover_level)
         {
            cell_state[neighbour] = LOWER_ADJACENT;
            lower_adjacent_cells.push(neighbour);
         }
         else
         {
            cell_state[neighbour] = UPPER_ADJACENT;
            upper_adjacent_cells.push(neighbour);
         }
      }
      else if (cell_state[neighbour] === COVERED)
      {
         if (covered_time[neighbour] > best_covered_time)
         {
            best_covered_time = covered_time[neighbour];
            last_covered_neighbour = neighbour;
            best_positions[cell*2] = position;
         }
         if (cell_elevation[neighbour] < lowest_elevation)
         {
            lowest_elevation = cell_elevation[neighbour];
            lowest_neighbour = neighbour;
            best_positions[cell*2+1] = position;
         }
      }
   }
   if (last_covered_neighbour != -1)
   {
      best_neighbours[cell*2] = last_covered_neighbour;
      in_degree[last_covered_neighbour]+=1;
   }
   if (lowest_neighbour != -1 && lowest_neighbour != last_covered_neighbour)
   {
      best_neighbours[cell*2+1] = lowest_neighbour;
      in_degree[lowest_neighbour]+=1;
   }
}

function fill_lakes()
{
   // Convenience helper used during initial setup.
   // Runs a drainage solve, then dumps lake water into sediment as a quick way
   // to pre-fill basins.
   initialize_drainage();
   calculate_drainage_graph();
   for (let cell = 0; cell<cell_count; cell++)
   {
      if (lake_index[cell] > 0)
      {
         sediment_thickness[cell] += lake_thickness[cell] + Math.random();
      }
   }
}

function run_drainage(dt, glaciate = false)
{
   // Route rainfall + upstream flow downstream, applying erosion/deposition.
   //
   // Execution model:
   // - The covering pass assigns each cell its downstream neighbours.
   // - `in_degree` counts how many upstream cells route into a given cell.
   // - We process cells in a topological order using `ready_to_drain_cells`.
   //
   // If you want to add a new step that depends on flow direction or discharge,
   // adding it inside this loop is usually the cleanest approach.
   let total_mass_eroded = 0;
   let total_mass_deposited = 0;
   let total_mass_off_edge = 0;
   for (let cell=0; cell<cell_count; cell++)
   {
      if (in_degree[cell] === 0) { ready_to_drain_cells.push(cell); }
   }
   while (ready_to_drain_cells.is_nonempty())
   {
      let cell = ready_to_drain_cells.pop();
      let best_neighbour = best_neighbours[cell*2];
      if (best_neighbour === -1) //we are at a sink
      {
         total_mass_off_edge += collected_sediment[cell];
         continue;
      }
      let best_position = best_positions[cell*2];

      let water = collected_water[cell];
      let sediment = collected_sediment[cell];
      let is_lake = (lake_thickness[cell] > lake_threshold);
      let elevation_here = cell_elevation[cell]-sea_level;
      let rainfall_here = Math.max(0.5+elevation_here/2000, 0);
      let cell_area = cell_size[cell] * pixel_length * pixel_length;
      water += rainfall_here * cell_area;

      let best_neighbour_level = cell_elevation[best_neighbour] + lake_thickness[best_neighbour];
      let best_distance = cell_neighbour_distance[best_position] * pixel_length;
      let level_here = cell_elevation[cell] + lake_thickness[cell];
      let best_gradient = (level_here - best_neighbour_level) / best_distance;
      if (best_gradient<0) best_gradient = 0; //?? I think to deal with shallow lakes
      /*if (!is_lake && level_here < best_neighbour_level)
      {
         console.log("lower neighbour");
      }*/

      // We will calculate erosion as though all water went to the best neighbour
      //let Q = Math.pow(water,2)/300000000;
      //let Q = Math.pow(water,1.7)/300000;
      //let Q = Math.pow(water, 0.8)*50;
      let Q = water;
      let Qm = Math.pow(Q, 0.5)*1000;

      //let under_ice = (glaciate && cell_elevation[cell]+lake_thickness[cell]-sea_level>300)
      //if (glaciate) Q = water;
      //if (!under_ice)

      let sediment_transport_capacity = (is_lake) ? 0 : fluvial_transport_coefficient * Q * best_gradient * dt;
      let equilibrium_slope = sediment / (fluvial_transport_coefficient * Q * dt);
      if (is_lake)
      {
         best_neighbour_level = cell_elevation[best_neighbour];
         equilibrium_slope = Math.min(equilibrium_slope, 0.5); //to stop sudden islands
      }
      let equilibrium_level = best_neighbour_level + equilibrium_slope * best_distance;
      // TODO - interpolate somehow between bedrock and sediment erosion
      if (sediment_transport_capacity > sediment && cell_elevation[cell] > equilibrium_level)
      {
         if (sediment_thickness[cell] > 0) Qm *= 2;
         //let erosion_distance = (sediment_thickness[cell]>0) ?
            //sediment_erosion_distance : bedrock_erosion_distance;
         //let mass_to_erode = (sediment_transport_capacity - sediment) *
            //   best_distance / erosion_distance;
         let mass_to_erode = Qm*best_gradient*dt * (sediment_transport_capacity - sediment)/sediment_transport_capacity;
         mass_to_erode *= erodibility[cell];
         let thickness_to_erode = mass_to_erode / (sediment_density * cell_area);
         if (cell_elevation[cell] - thickness_to_erode < equilibrium_level)
         {
            thickness_to_erode = cell_elevation[cell] - equilibrium_level;
            mass_to_erode = thickness_to_erode * cell_area * sediment_density;
         }
         sediment_thickness[cell] -= thickness_to_erode;
         if (sediment_thickness[cell]<0)
         {
            bedrock_thickness[cell] += sediment_thickness[cell];
            sediment_thickness[cell] = 0;
         }
         sediment += mass_to_erode;
         total_mass_eroded += mass_to_erode;
      }
      else if (sediment_transport_capacity < sediment && cell_elevation[cell] < equilibrium_level)
      {
         let mass_to_deposit = (sediment - sediment_transport_capacity) *
               best_distance / sedimentation_distance;
         //if (under_ice) { mass_to_deposit *= 0.5; }
         let thickness_to_deposit = mass_to_deposit / (sediment_density * cell_area);
         if (cell_elevation[cell] + thickness_to_deposit > equilibrium_level)
         {
            thickness_to_deposit = equilibrium_level - cell_elevation[cell];
            mass_to_deposit = thickness_to_deposit * cell_area * sediment_density;
         }
         sediment_thickness[cell] += thickness_to_deposit;
         sediment -= mass_to_deposit;
         total_mass_deposited += mass_to_deposit;
      }

      /*if (under_ice)
      {
         let glacier_erosion_factor = 0.00000001;
         let thickness_to_erode = glacier_erosion_factor * water;
         let mass_to_erode = thickness_to_erode * sediment_density * cell_area;
         sediment_thickness[cell] -= thickness_to_erode;
         if (sediment_thickness[cell]<0)
         {
            bedrock_thickness[cell] += sediment_thickness[cell];
            sediment_thickness[cell] = 0;
         }
         sediment += mass_to_erode;
      }*/

      in_degree[best_neighbour] -= 1;
      if (in_degree[best_neighbour] === 0) { ready_to_drain_cells.push(best_neighbour); }

      let second_neighbour = best_neighbours[cell*2+1];
      if (second_neighbour === -1 || second_neighbour === best_neighbour) //ie. only one neighbour
      {
         discharge[best_position] = water;
         collected_water[best_neighbour] += water;
         collected_sediment[best_neighbour] += sediment;
      }
      else
      {
         in_degree[second_neighbour] -= 1;
         if (in_degree[second_neighbour] === 0) { ready_to_drain_cells.push(second_neighbour); }
         let second_position = best_positions[cell*2+1];

         let first_fraction = 1.0;//0.667;
         let second_fraction = 0.0;//0.333;

         if (!is_lake)
         {
            let second_neighbour_level = cell_elevation[second_neighbour] + lake_thickness[second_neighbour];
            let second_distance = cell_neighbour_distance[second_position] * pixel_length;
            let second_gradient = (level_here - second_neighbour_level) / second_distance;// (multiplier * pixel_length);
            if (second_gradient > best_gradient)
            {
               second_gradient = best_gradient;
            }
            best_gradient += 0.000000000001;
            second_gradient = Math.pow(second_gradient/best_gradient, 10);
            let total_gradient = 1 + second_gradient; //assume we have divided best_gradient by itself
            first_fraction = 1/total_gradient;
            second_fraction = second_gradient/total_gradient;
            if (water < 4*km3) // 2 times Thames?
            {
               first_fraction = 1.0;
               second_fraction = 0.0;
            }
         }

         discharge[best_position] = water * first_fraction;
         collected_water[best_neighbour] += water * first_fraction;
         collected_sediment[best_neighbour] += sediment * first_fraction;

         discharge[second_position] = water * second_fraction;
         collected_water[second_neighbour] += water * second_fraction;
         collected_sediment[second_neighbour] += sediment * second_fraction;
      }


      //if (cell_elevation[cell] > sea_level+1500)
      {

         let delta = Math.max(0, cell_elevation[cell] - cell_elevation[best_neighbour]-2*pixel_length);
         bedrock_thickness[cell] -= delta/2;
         sediment_thickness[best_neighbour] += delta/2;

         let k = 0*0.5* dt / 1000000;
         delta = (cell_elevation[cell] - cell_elevation[best_neighbour])*k;
         bedrock_thickness[cell] -= delta;
         sediment_thickness[best_neighbour] += delta;
      }

      //if ((lake_thickness[best_neighbour]>1  && lake_thickness[best_neighbour]<120)
         //&& (cell_elevation[cell]>sea_level-100 && cell_elevation[cell] < sea_level+100))
      if (lake_index[best_neighbour] === 0)
      {
         let k = landscape_diffusion_coefficient * dt / 1000000;
         k *= 3;
         let delta = (cell_elevation[cell] - cell_elevation[best_neighbour])*k;
         if (sediment_thickness[cell] < delta) delta = sediment_thickness[cell];
         sediment_thickness[cell] -= delta;
         sediment_thickness[best_neighbour] += delta;
      }
   }
   /*low_readout.innerHTML = "eroded: "+total_mass_eroded.toExponential(2) +
         "  &nbsp  deposited / off edge: " +
        (total_mass_deposited+total_mass_off_edge).toExponential(2) + " = " +
        total_mass_deposited.toExponential(2) + " / " +
        total_mass_off_edge.toExponential(2);*/

}
