/* Copyright (C) 2011,2012 Chi-kwan Chan
   Copyright (C) 2011,2012 NORDITA

   This file is part of fg2.

   Fg2 is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Fg2 is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with fg2.  If not, see <http://www.gnu.org/licenses/>. */

#include "fg2.h"
#include "deriv.h"
#include <eqns.cuh>

static __device__ void copy(R *dst, const R *src, const Z h, const Z s)
{
  const Z w = (blockDim.x + ORDER) * NVAR;
  const Z n =  blockDim.x * blockDim.y;
  const Z m = (w - 1) / n + 1;
  const Z l = threadIdx.y * blockDim.x + threadIdx.x;

  for(Z i = 0; i < h; ++i) {
    for(Z j = 0; j < m; ++j) {
      const Z k = j * n + l;
      if(k < w) dst[i * w + k] = src[i * s + k];
    }
    __syncthreads();
  }
}

static __device__ void ssum(R *dst, const R *src, const R beta, const Z s)
{
  const Z w =  blockDim.x * NVAR;
  const Z n =  blockDim.x * blockDim.y;
  const Z m = (w - 1) / n + 1;
  const Z l = threadIdx.y * blockDim.x + threadIdx.x;

  for(Z i = 0; i < blockDim.y; ++i) {
    for(Z j = 0; j < m; ++j) {
      const Z k = j * n + l;
      if(k < w) dst[i * s + k] = beta * dst[i * s + k] + src[i * w + k];
    }
    __syncthreads();
  }
}

static __global__ void kernel(R *out, const R *in, const R t, const R beta,
                              Z n1, Z n2, const Z s)
{
  extern __shared__ R shared[]; // dynamic shared variable

  Z i1 = threadIdx.y;
  Z i2 = threadIdx.x;
  {
    const Z m1 = blockDim.y * ((n1 - 1) / (gridDim.y * blockDim.y) + 1);
    const Z m2 = blockDim.x;

    const Z offset = blockIdx.y * m1 * s + blockIdx.x * m2 * NVAR;
    out += offset;
    in  += offset + HALF * (s - NVAR);

    n1 = (blockIdx.y == gridDim.y - 1) ? n1 - blockIdx.y * m1 : m1;
    n2 = (blockIdx.x == gridDim.x - 1) ? n2 - blockIdx.x * m2 : m2;

    i1 += blockIdx.y * m1;
    i2 += blockIdx.x * m2;
  }
  const Z w = (blockDim.x + ORDER) * NVAR;
  const Z h = (n1 - 1) / blockDim.y + 1;
  const Z l = threadIdx.y * blockDim.x + threadIdx.x;

  const S *active = (S *)(shared + (HALF + threadIdx.y) * w
                                 + (HALF + threadIdx.x) * NVAR);

  // Modified rolling cache (Micikevicius 2009)
  for(Z i = 0; i < h; ++i) {

    if(i == 0) copy(shared, in     - ORDER      * s, ORDER, s);
    else       copy(shared, shared + blockDim.y * w, ORDER, w);

    copy(shared + ORDER * w, in + i * blockDim.y * s, blockDim.y, s);

    const S f = eqns(active, i1 + i * blockDim.y, i2, blockDim.x + ORDER);
    __syncthreads();

    ((S *)shared)[l] = f;
    __syncthreads();

    ssum(out + i * blockDim.y * s, shared, beta, s);
  }
}

void kick(R *v, const R *x, const R t, const R b)
{
  using namespace global;
  const dim3 Gsz(g2, g1);
  const dim3 Bsz(b2, b1);
  kernel<<<Gsz, Bsz, sz>>>(v, x, t, b, n1, n2, s);
}

// CUDA cannot really understand constant device variables across many
// files.  Therefore we are forced to put the Euler kernel here.

#include <Euler.cuh>

static __global__ void kernel(R *u, const R dt, const Z n, const Z s)
{
  const Z i = blockIdx.y;
  const Z j = blockIdx.x * blockDim.x + threadIdx.x;
  if(j < n) first_order((S *)(u + i * s + j * NVAR), dt, i, j);
}

void Euler(R *x, const R dt)
{
  using namespace global;
  const dim3 Gsz((n2 - 1) / BSZ + 1, n1);

#ifdef RANDOMIZE
  RANDOMIZE
#endif

  kernel<<<Gsz, BSZ>>>(x, dt, n2, s);
}
