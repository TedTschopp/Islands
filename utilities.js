"use strict";

/*
Utilities

This file provides a small grab-bag of helpers used across the simulation.

- Morton/Z-order helpers:
  - `interleave_bits` / `deinterleave_bits` are used when `morton_order` is enabled
    in `create_voronoi_tiling` (see `voronoi.js`). In this project morton order is
    optional and usually disabled.

- Noise:
  - `SNoise` wraps the global `Noise.simplex2` (from `noise.js`) to provide:
    - anisotropic wavelength (different wavelength in s/t)
    - rotation
    - multiple octaves with controlled rolloff

- Small data structures used by drainage (`drainage.js`):
  - `Stack` (LIFO)
  - `RandomQueue` (biased random pop from the "old" end; encourages organic growth)
  - `PriorityQueue` (binary heap keyed by an external `values[]` array)

Implementation notes:
- This code intentionally uses TypedArrays for performance.
- `blend_colours` appears twice (legacy duplication). The second definition wins.
*/

function interleave_bits(str, k) //str is a 2k-bit string
{
   let left = str>>k;
   let right = str; //will not use leftmost bits
   let output = 0;
   for (let i=0; i<k; i++)
   {
      output = set_bit(output, 2*i, left & 1);
      output = set_bit(output, 2*i+1, right & 1);
      left = left>>1;
      right = right>>1;
   }
   return output;
}

function deinterleave_bits(str, k)
{
   let output = 0;
   for (let i=0; i<k; i++)
   {
      let left_bit = str & 1;
      str = str >> 1;
      let right_bit = str & 1;
      str = str >> 1;
      output = set_bit(output, k+i, left_bit);
      output = set_bit(output, i, right_bit);
   }
   return output;
}

function set_bit(str, k, b) { return str | (b<<k); }

function random_integer(n)
{
   return Math.floor(Math.random() * n);
}

var output = [0,0,0];
function blend_colours(a, b, s)
{
  for (var i=0; i<3; i++) output[i] = a[i]*(1-s) + b[i]*s;
  return output;
}

class SNoise
{
  constructor(wavelength, cross_wavelength=-1, angle=0)
  {
    this.wavelength_t = wavelength;
    if (cross_wavelength == -1) this.wavelength_s = wavelength;
    else this.wavelength_s = cross_wavelength;
    this.offset_s = Math.random()*100;
    this.offset_t = Math.random()*100;
    this.theta = 2*Math.PI*angle/360;
    this.sin_theta = Math.sin(this.theta);
    this.cos_theta = Math.cos(this.theta);
    this.octaves_count = 1;
    this.wavenumber_rolloff = 2;
    this.amplitude_rolloff = 0.5;
    this.rms = 1;
  }

  octaves(count, amplitude_rolloff = 0.5, wavenumber_rolloff = 2)
  {
    this.octaves_count = count;
    this.wavenumber_rolloff = wavenumber_rolloff;
    this.amplitude_rolloff = amplitude_rolloff;
    var square = 0;
    for (var i=0; i<count; i++) square += Math.pow(amplitude_rolloff, 2*i);
    this.rms = Math.sqrt(square);
  }

  get(x, y)
  {
    var unscaled_s = x * this.cos_theta + y * this.sin_theta;
    var unscaled_t = x * this.sin_theta - y * this.cos_theta;
    var t = unscaled_t / this.wavelength_t;
    var s;
    if (this.wavelength_s == 0) s = 0;
    else s = unscaled_s / this.wavelength_s;
    //return Noise.simplex2(s+this.offset_s, t+this.offset_t)
    return this.get_basic_noise(s,t);
  }

  get_basic_noise(s,t)
  {
    var wavenumber = 1;
    var amplitude = 1;
    var sum = 0;
    for (var i=0; i<this.octaves_count; i++)
    {
      var signal = Noise.simplex2(s*wavenumber+(i+1)*this.offset_s, t*wavenumber+(i+1)*this.offset_t);
      sum += signal * amplitude;
      amplitude *= this.amplitude_rolloff;
      wavenumber *= this.wavenumber_rolloff;
    }
    return sum / this.rms;
  }
}

var output = [0,0,0];
function blend_colours(a, b, s)
{
  for (var i=0; i<3; i++) output[i] = a[i]*(1-s) + b[i]*s;
  return output;
}

class Stack
{
  constructor(n)
  {
    this.indices = new Uint32Array(n);
    this.pointer = -1;
  }

  clear()
  {
    this.pointer = -1;
  }

  is_nonempty()
  {
    return (this.pointer>-1);
  }

  push(index)
  {
    this.pointer +=1;
    this.indices[this.pointer] = index;
  }

  pop()
  {
    let output = this.indices[this.pointer];
    this.pointer -= 1;
    return output;
  }
}

class RandomQueue //returns randomly one of the oldest elements in the queue
{
   constructor(n)
   {
      this.indices = new Uint32Array(n);
      this.start_pointer = 0;
      this.end_pointer = 0;
   }

   clear()
   {
      this.start_pointer = 0;
      this.end_pointer = 0;
   }

   is_nonempty()
   {
      return (this.end_pointer > this.start_pointer);
   }

   push(index)
   {
      //if (this.pointer == n) { this.compact(); }
      this.indices[this.end_pointer]=index;
      this.end_pointer += 1;
   }

   pop()
   {
      let offset = Math.floor(Math.random() * 0.1 * (this.end_pointer - this.start_pointer));
      let position = this.start_pointer + offset;
      if (position >= this.end_pointer) position = this.start_pointer;
      let output = this.indices[position];
      this.indices[position] = this.indices[this.start_pointer];
      this.start_pointer += 1;
      return output;
   }
}

class PriorityQueue // TODO - better to store values in heap, next to index?
{
  constructor(n, values)
  {
    this.indices_heap = new Uint32Array(n+1);
    this.values = values;
    this.values_heap = new Float32Array(n+1);
    this.heap_pointer = 0;
    this.max_heap_size = 0;
  }

  better(x, y)  // the best thing goes to the top of the heap
  {
    return (x<y);
  }

  clear()
  {
    //console.log(this.max_heap_size);
    this.max_heap_size = 0;
    this.heap_pointer = 0;
  }

  is_nonempty()
  {
    return (this.heap_pointer>0);
  }

  push(index)
  {
    this.heap_pointer += 1;
    if (this.heap_pointer > this.max_heap_size) this.max_heap_size = this.heap_pointer;
    let value = this.values[index];
    let pos = this.heap_pointer;
    let parent;
    let parent_value;
    while (pos>1)
    {
      parent = this.indices_heap[pos>>1];
      parent_value = this.values_heap[pos>>1];
      if (this.better(value, parent_value))
      {
        this.indices_heap[pos] = parent;
        this.values_heap[pos] = parent_value;
        pos = pos>>1;
      }
      else break;
    }
    this.indices_heap[pos] = index;
    this.values_heap[pos] = value;
  }

  pop()
  {
    let output = this.indices_heap[1];
    let pos = 1;
    let left_child, right_child, left_value, right_value;
    while (true) // propagate hole down to leaf
    {
      let left = pos<<1;
      if (this.heap_pointer < left) //pos has no children
      {
        break;
      }
      else if (this.heap_pointer == left) //only left child
      {
        left_child = this.indices_heap[left];
        left_value = this.values_heap[left];
        this.indices_heap[pos] = left_child; //child is smaller
        this.values_heap[pos] = left_value;
        pos = left;
        continue;
      }
      else //two children
      {
        left_child = this.indices_heap[left];
        left_value = this.values_heap[left];
        right_child = this.indices_heap[left+1];
        right_value = this.values_heap[left+1];
        if (this.better(left_value, right_value))
        {
          this.indices_heap[pos] = left_child; //swap with left child
          this.values_heap[pos] = left_value;
          pos = left;
          continue;
        }
        else
        {
          this.indices_heap[pos] = right_child; //swap with right child
          this.values_heap[pos] = right_value;
          pos = left+1;
          continue;
        }
      }
    } // hole is now at a leaf, pos

    // now move last entry into hole, and propagate up
    let index = this.indices_heap[this.heap_pointer];
    let value = this.values_heap[this.heap_pointer];
    this.heap_pointer -= 1;
    let parent;
    let parent_value;
    while (pos>1)
    {
      parent = this.indices_heap[pos>>1];
      parent_value = this.values_heap[pos>>1];
      if (this.better(value, parent_value))
      {
        this.indices_heap[pos] = parent;
        this.values_heap[pos] = parent_value;
        pos = pos>>1;
      }
      else break;
    }
    this.indices_heap[pos] = index;
    this.values_heap[pos] = value;

    return output;
  }
}
