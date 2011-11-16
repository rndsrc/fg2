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

static __global__ void kernel(R *x, const R *v, const R dt,
                              const Z n, const Z s)
{
  const Z j = blockIdx.x * blockDim.x + threadIdx.x;
  if(j < n) {
    const Z h = blockIdx.y * s + j;
    x[h] += dt * v[h];
  }
}

void drift(R *x, const R *v, const R dt)
{
  const Z n = NVAR * global::n2;
  const dim3 Gsz((n - 1) / BSZ + 1, global::n1);
  kernel<<<Gsz, BSZ>>>(x, v, dt, n, global::s);
}
