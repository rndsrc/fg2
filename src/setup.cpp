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

#include <cstdlib>
#include <limits>
#include <cuda_runtime.h>
#include "fg2.h"

namespace global {
  Z n1, n2, s;
  R *u, *v, *host = NULL;
  Z g1, g2, b1, b2, sz;
  double flops, bps;
}

static void done(void)
{
  if(global::host)
   free(global::host - HALF * (global::s + NVAR));
  cudaFree(global::v - HALF * (global::s + NVAR));
  cudaFree(global::u - HALF * (global::s + NVAR));
}

int setup(Z n1, Z n2)
{
  if(atexit(done)) abort();

  Z m1 = (global::n1 = n1) + ORDER;
  Z m2 = (global::n2 = n2) + ORDER;

  // Grid and block sizes for rolling cache kernel
  int d; cudaGetDevice(&d);
  cudaDeviceProp dev; cudaGetDeviceProperties(&dev, d);
  for(Z h = 0, i = 0, j = 1; i < j && i++ * j < BSZ; ) {
    j = dev.sharedMemPerBlock / (sizeof(S) * (ORDER + i)) - ORDER;
    if(j > dev.maxThreadsPerBlock) j = dev.maxThreadsPerBlock;
    else j = (j / dev.warpSize) * dev.warpSize;
    if(i * j > h) {
      h = i * j;
      global::b1 = i;
      global::b2 = j;
      global::sz = sizeof(S) * (ORDER + i) * (ORDER + j);
    }
  }
  global::g2 = (n2 - 1) / global::b2 + 1;
  global::g1 = (dev.multiProcessorCount - 1) / global::g2 + 1;

  // State variable
  void *u; size_t upitch;
  if(cudaSuccess != cudaMallocPitch(&u, &upitch, sizeof(S) * m2, m1) ||
     upitch % sizeof(R)) return 0;
  global::s = upitch / sizeof(R);
  global::u = (R *)u + HALF * (global::s + NVAR);

  // Storage for finite difference or swap space for finite volume
  void *v; size_t vpitch;
  if(cudaSuccess != cudaMallocPitch(&v, &vpitch, sizeof(S) * m2, m1) ||
     vpitch % sizeof(R)) return 0;
  if(vpitch != upitch) return 0;
  global::v = (R *)v + HALF * (global::s + NVAR);

  // Allocate host memory
  Z  n = global::s * m1;
  R *h = (R *)malloc(sizeof(R) * n);
  if(NULL == h) return 0;
  global::host = h + HALF * (global::s + NVAR);

  // Initialize all arrays to -FLT_MAX or -DBL_MAX
  for(Z i = 0; i < n; ++i) h[i] = 0.0; // TODO: set boundary conditions
                               // -std::numeric_limits<R>::max();
  cudaMemcpy(u, h, sizeof(R) * n, cudaMemcpyHostToDevice);
  cudaMemcpy(v, h, sizeof(R) * n, cudaMemcpyHostToDevice);

  // Compute floating point operation and bandwidth per step
  global::flops = 3 * ((n1 * n2) * (181 + NVAR * 2.0)); // assume FMA
  global::bps   = 3 * ((m1 * m2) *         sizeof(S) * 8 * 1.0 +
                       (n1 * n2) *         sizeof(S) * 8 * 5.0 +
                       (m1 + m2) * ORDER * sizeof(S) * 8 * 2.0);

  // Return size of device memory
  return 2 * sizeof(R) * n;
}
