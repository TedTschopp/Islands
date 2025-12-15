"use strict"

/*
Voronoi tiling + triangle mesh construction

This file constructs the simulation's spatial discretization.

The core entrypoint is:

   create_voronoi_tiling(cellWidth, cellHeight, multiplier, mortonBits=0, relaxationCount=2)

Conceptually:
- You choose a regular grid of logical "cells" (`cellWidth` x `cellHeight`).
- Each cell is assigned a random site position inside its pixel tile.
- Every pixel in the high-res raster (`pixelWidth` x `pixelHeight`) is assigned to
   its nearest site (`pixels[p] = cell`).
- The sites are relaxed a few times (Lloyd-like) to reduce clustering.
- Adjacency is recovered by scanning pixel edges and recording cell transitions.
- Triangles are built by connecting each cell to consecutive neighbours around it.

Important globals produced (used by simulation + rendering):

- Dimensions:
   - `cell_width`, `cell_height`, `cell_count`
   - `pixel_width = cell_width * multiplier`, `pixel_height = cell_height * multiplier`

- Per-cell geometry:
   - `cell_x[cell]`, `cell_y[cell]`: site location in pixel coordinates
   - `cell_size[cell]`: number of pixels owned by the cell (approx area)

- Neighbourhood:
   - `max_neighbours`: fixed upper bound for neighbour list storage
   - `cell_neighbour_count[cell]`: number of neighbours for this cell
   - `cell_neighbours[cell*max_neighbours + i]`: neighbour cell id
   - `cell_neighbour_distance[...]`: Euclidean distance between sites (in pixel units)
   - `cell_neighbours_coded[...]`: compact neighbour offsets (used by `get_neighbour`)

- Triangle mesh for WebGL:
   - `triangle_count`
   - `triangle_cells[tri*3 + j]`: the 3 cell ids that define triangle vertices
   - `triangle_vertices[tri*6 + 2*j + {0,1}]`: the (x,y) coords for each vertex

If you want a different resolution:
- Change the call in `index.html`.
- Be mindful that `pixel_width`/`pixel_height` (and thus the canvas) scale with
   `multiplier`.
*/

let max_neighbours = 12;
let pixel_width, pixel_height, pixel_count;
let multiplier;
let cell_count, cell_width, cell_height;
let cell_x, cell_y, cell_size;
let cell_neighbour_count, cell_neighbours, cell_neighbour_distance;
let triangle_cells, triangle_vertices, triangle_count;
let get_nearest_cell;
let morton_order, morton_bits;
let cell_neighbours_coded;
let cell_colour, pixels;

function create_voronoi_tiling(width_in, height_in, multiplier_in, morton_bits_in = 0, relaxation_count = 2)
{
   let cell_total_x, cell_total_y;
   let neighbour_counts = new Uint32Array(max_neighbours);
   if (morton_bits_in != 0)
   {
      morton_bits = morton_bits_in;
      morton_order = true;
   }

   function create(width_in, height_in, multiplier_in)
   {
      multiplier = multiplier_in;
      cell_width = width_in;
      cell_height = height_in;
      pixel_width = width_in * multiplier;
      pixel_height = height_in * multiplier;
      pixel_count = pixel_width * pixel_height;
      cell_count = width_in * height_in;

      pixels = new Uint32Array(pixel_count);
      cell_x = new Float32Array(cell_count);
      cell_y = new Float32Array(cell_count);
      cell_total_x = new Float32Array(cell_count);
      cell_total_y = new Float32Array(cell_count);
      cell_size = new Uint16Array(cell_count);
      cell_neighbour_count = new Uint8Array(cell_count);
      cell_neighbours = new Uint32Array(cell_count*max_neighbours);
      cell_neighbours_coded = new Uint8Array(cell_count * max_neighbours);
      cell_neighbour_distance = new Float32Array(cell_count*max_neighbours);
      cell_colour = new Uint32Array(cell_count*3);
      //row_to_morton = new Uint32Array(cell_count);
      //morton_to_row = new Uint32Array(cell_count);


      for (let i=0; i<pixel_count; i++) { pixels[i] = -1; }

      for (let j=0; j<height_in; j++)
      {
         for (let i=0; i<width_in; i++)
         {
            let cell = j * width_in + i;
            if (morton_order) { cell = interleave_bits(cell, morton_bits); }
            let y = j * multiplier + Math.random()*multiplier;
            let x = i * multiplier + Math.random()*multiplier;
            cell_x[cell] = x;
            cell_y[cell] = y;
         }
      }

      for (let i=0; i<cell_count*3; i++)
      {
         cell_colour[i] = 50+random_integer(200);
      }
   }

   function tile()
   {
      for (let cell=0; cell<cell_count; cell++)
      {
         cell_total_x[cell] = cell_total_y[cell] = cell_size[cell] = 0;
      }

      for (let y=0; y<pixel_height; y++)
      {
         for (let x=0; x<pixel_width; x++)
         {
            let cell = get_nearest_cell(x,y);
            pixels[y * pixel_width + x] = cell;
            cell_total_x[cell] += x;
            cell_total_y[cell] += y;
            cell_size[cell] += 1;
         }
      }
   }

   function relax(factor = 1.0)
   {
      for (let cell=0; cell<cell_count; cell++)
      {
         let centroid_x = cell_total_x[cell]/cell_size[cell];
         let centroid_y = cell_total_y[cell]/cell_size[cell];

         cell_x[cell] = (1-factor)*cell_x[cell] + factor*centroid_x;
         cell_y[cell] = (1-factor)*cell_y[cell] + factor*centroid_y;
      }
   }

   get_nearest_cell = function(x,y)
   {
      const search_distance = 1;
      const centre_i = Math.floor(x/multiplier);
      const centre_j = Math.floor(y/multiplier);
      let smallest_distance_squared = pixel_count;
      let nearest_cell = -1;
      for (let dj = -search_distance; dj<=search_distance; dj++)
      {
         for (let di = -search_distance; di<=search_distance; di++)
         {
            let i = centre_i + di;
            let j = centre_j + dj;
            let cell = j * cell_width + i;
            if (cell < 0 || cell >= cell_count) continue;
            if (morton_order) { cell = interleave_bits(cell, morton_bits); }
            let delta_x = cell_x[cell] - x;
            let delta_y = cell_y[cell] - y;
            let distance_squared = delta_x*delta_x + delta_y*delta_y;
            if (distance_squared < smallest_distance_squared)
            {
               smallest_distance_squared = distance_squared;
               nearest_cell = cell;
            }
         }
      }
      return nearest_cell;
   }

   function find_neighbours()
   {
      for (let i=0; i<cell_count; i++) cell_neighbour_count[i] = 0;
      for (let i=0; i<cell_count*max_neighbours; i++) cell_neighbours[i] = 0;

      for (let y=0; y<pixel_height; y++)
      {
         for (let x=0; x<pixel_width; x++)
         {
            let i = y * pixel_width + x;
            let this_cell = pixels[i];
            if (y<pixel_height-1)
            {
               let down_cell = pixels[i+pixel_width];
               if (down_cell != this_cell
                  && (cell_neighbour_count[this_cell] < max_neighbours
                  && cell_neighbour_count[down_cell] < max_neighbours))
               {
                  add_neighbour(this_cell, down_cell);
                  add_neighbour(down_cell, this_cell);
               }
            }
            if (x<pixel_width-1)
            {
               let right_cell = pixels[i+1];
               if (right_cell != this_cell
                  && (cell_neighbour_count[this_cell] < max_neighbours
                  && cell_neighbour_count[right_cell] < max_neighbours))
               {
                  add_neighbour(this_cell, right_cell);
                  add_neighbour(right_cell, this_cell);
               }
            }
            if (y<pixel_height-1 && x<pixel_width-1)
            {
               let down_right_cell = pixels[i+pixel_width+1];
               if (down_right_cell != this_cell
                  && (cell_neighbour_count[this_cell] < max_neighbours
                  && cell_neighbour_count[down_right_cell] < max_neighbours))
               {
                  add_neighbour(this_cell, down_right_cell);
                  add_neighbour(down_right_cell, this_cell);
               }
            }
         }
      }

      process_neighbours();
   }

   function add_neighbour(this_cell, other_cell)
   {
      let neighbour_count = cell_neighbour_count[this_cell];
      let position = this_cell * max_neighbours;
      for (let i = position; i<position+neighbour_count; i++)
      {
         if (cell_neighbours[i] === other_cell) return;
      }
      cell_neighbours[position+neighbour_count] = other_cell;
      cell_neighbour_count[this_cell] +=1;
   }

   function process_neighbours()
   {
      for (let cell=0; cell<cell_count; cell++)
      {
         neighbour_counts[cell_neighbour_count[cell]] += 1;
         let indices = [];
         let neighbours = [];
         let distances = [];
         let angles = [];
         for (let i=0; i<cell_neighbour_count[cell]; i++)
         {
            let position = cell*max_neighbours + i;
            let neighbour = cell_neighbours[position];
            let dx = cell_x[neighbour] - cell_x[cell];
            let dy = cell_y[neighbour] - cell_y[cell];
            let distance = Math.sqrt(dx*dx+dy*dy);
            let angle = Math.atan2(dy, dx);
            indices.push(i);
            neighbours.push(neighbour);
            distances.push(distance);
            angles.push(angle);
         }
         indices.sort(function(i,j) { return (angles[i]>angles[j]);} );
         for (let i=0; i<cell_neighbour_count[cell]; i++)
         {
            let position = cell*max_neighbours + i;
            let index = indices[i];
            cell_neighbours[position] = neighbours[index];
            cell_neighbour_distance[position] = distances[index];
         }
      }
      /*for (let i=0; i<max_neighbours; i++)
      {
         console.log(i+": "+neighbour_counts[i]);
      }*/
   }

   function create_triangles() //TODO remove edge triangles (detect by large angles?)
   {
      let hash_set = {};
      let triangle_cells_list = [];
      triangle_count = 0;
      let duplicate_triangle_count = 0;
      for (let cell=0; cell<cell_count; cell++)
      {
         for (let i=0; i<cell_neighbour_count[cell]; i++)
         {
            let position1 = cell*max_neighbours + i;
            let position2 = cell*max_neighbours + (i+1)%cell_neighbour_count[cell];
            let neighbour1 = cell_neighbours[position1];
            let neighbour2 = cell_neighbours[position2];
            //let sorted = [cell, neighbour1, neighbour2].sort();
            let a,b,c,d,e;
            if (cell<neighbour1) { d=cell; e=neighbour1; }
            else { d = neighbour1; e = cell; }
            if (neighbour2 < d) { a = neighbour2; b=d; c=e; }
            else if (neighbour2 < e) { a = d; b = neighbour2; c = e; }
            else { a=d; b=e; c=neighbour2; }
            let sorted = a+" "+b+" "+c;
            if (sorted in hash_set)
            {
               duplicate_triangle_count += 1;
               continue;
            }
            hash_set[sorted] = true;
            triangle_cells_list.push(cell, neighbour1, neighbour2);
            triangle_count += 1;
         }
      }
      //console.log(triangle_count, duplicate_triangle_count);
      triangle_cells = new Uint32Array(triangle_cells_list);
      triangle_vertices = new Float32Array(triangle_count*6);
      for (let i=0; i<triangle_count*3; i++)
      {
         let cell = triangle_cells[i];
         triangle_vertices[2*i] = cell_x[cell];
         triangle_vertices[2*i+1] = cell_y[cell];
      }
   }

   function get_average_neighbour_count()
   {
      let total = 0;
      for (let cell=0; cell<cell_count; cell++) total += cell_neighbour_count[cell];
      return total/cell_count;
   }

   function setup_coded_neighbours()
   {
      for (let cell=0; cell<cell_count; cell++)
      {
         for (let i=0; i<cell_neighbour_count[cell]; i++)
         {
            let position = cell*max_neighbours+i;
            cell_neighbours_coded[position] = get_neighbour_code(cell, cell_neighbours[position]);
         }
      }
   }

   create(width_in, height_in, multiplier_in);
   for (let i=0; i<relaxation_count; i++)
   {
      tile();
      relax(1.1);
   }
   tile();
   find_neighbours();
   setup_coded_neighbours();
   create_triangles();
}

function get_neighbour(cell, code)
{
   if (morton_order) { cell = deinterleave_bits(cell, morton_bits); }
   let x = cell & 255;
   let y = cell >> 8;
   let i = code & 15; // bottom 4 bits
   let j = code >> 4; // top 4 bits
   i = i - 8;
   j = j - 8;
   let neighbour_x = x+i;
   let neighbour_y = y+j;
   let out_cell = neighbour_y * cell_width + neighbour_x;
   if (morton_order) { return interleave_bits(out_cell, morton_bits); }
   else { return out_cell; }
}

function get_neighbour_code(cell, neighbour)
{
   if (morton_order)
   {
      cell = deinterleave_bits(cell, morton_bits);
      neighbour = deinterleave_bits(neighbour, morton_bits);
   }
   let x = cell & 255;
   let y = cell >> 8;
   let neighbour_x = neighbour & 255;
   let neighbour_y = neighbour >> 8;
   let i = (neighbour_x - x)+8;
   let j = (neighbour_y - y)+8;
   //if (Math.abs(i-8)>5 || Math.abs(j-8)>5) console.log("NEIGHBOUR TOO FAR");
   return (j<<4)+i;
}
