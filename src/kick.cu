/* Copyright (C) 2011 Chi-kwan Chan
   Copyright (C) 2011 NORDITA

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
#include STRINGIZE(EQNS)

static __device__ void copy(R *dst, const R *src, const Z n)
{
  const Z m = (n - 1) / blockDim.x + 1;
  for(Z i = 0; i < m; ++i) {
    const Z j = i * blockDim.x + threadIdx.x;
    if(j < n) dst[j] = src[j];
  }
}

static __device__ void ssum(R *dst, const R *src, const R beta, const Z n)
{
  const Z m = (n - 1) / blockDim.x + 1;
  for(Z i = 0; i < m; ++i) {
    const Z j = i * blockDim.x + threadIdx.x;
    if(j < n) dst[j] = beta * dst[j] + src[j];
  }
}

static __global__ void kernel(R *v, const R *x, const R t, const R beta,
                              Z n1, Z n2, const Z s)
{
  extern __shared__ R shared[]; // dynamic shared variable
  {
    const Z n = blockDim.x;               // inner size
    const Z m = (n1 - 1) / gridDim.y + 1; // outer size
    const Z l = blockIdx.y * m * s + blockIdx.x * n * NVAR;
    // Offset global arrays to the correct locations
    v += l;
    x += l + HALF * (s - NVAR);
    // Redefine local size
    n1 = (blockIdx.y == gridDim.y - 1) ? n1 - blockIdx.y * m : m;
    n2 = (blockIdx.x == gridDim.x - 1) ? n2 - blockIdx.x * n : n;
  }
  const Z Count = (n2 + ORDER) * NVAR;
  const Z count = (n2        ) * NVAR;

  R *in     = shared + (threadIdx.y + ORDER) * Count;
  R *active = shared + (threadIdx.y + HALF ) * Count + HALF * NVAR;
  R *out    = shared + (threadIdx.y        ) * count;

  // Modified rolling cache (Micikevicius 2009)
  for(Z i = threadIdx.y; i < n1; i += blockDim.y) {
    if(i < blockDim.y) // pre-load cache
      for(Z j = threadIdx.y; j < ORDER; j += blockDim.y)
        copy(shared + j * Count, x + (j - ORDER) * s, Count);
    else // shift cache
      for(Z j = threadIdx.y; j < ORDER; j += blockDim.y)
        copy(shared + j * Count, shared + (j + blockDim.y) * Count, Count);
    copy(in, x + i * s, Count);
    __syncthreads();

    const Z i1 = blockIdx.y * n1 + i;
    const Z i2 = blockIdx.x * blockDim.x + threadIdx.x;
    const S f  = eqns((S *)active + threadIdx.x, i1, i2, n2 + ORDER);
    __syncthreads();

    if(threadIdx.x < n2) *((S *)out + threadIdx.x) = f;
    __syncthreads();

    ssum(v + i * s, out, beta, count);
    __syncthreads();
  }
}

void kick(R *v, const R *x, const R t, const R b)
{
  using namespace global;
  const dim3 Gsz(g2, g1);
  const dim3 Bsz(b2, b1);
  kernel<<<Gsz, Bsz, sz>>>(v, x, t, b, n1, n2, s);
}
