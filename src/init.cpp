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

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cuda_runtime.h>
#include "fg2.h"
#include <init.h>

void init(const char *name)
{
  S (*func)(R, R) = pick(name);

  using namespace global;
  for(Z i = 0; i < n1; ++i) {
    const R x = l1 * (i + 0.5) / n1;
    for(Z j = 0; j < n2; ++j) {
      const R y = l2 * (j + 0.5) / n2;
      ((S *)host)[i * n2 + j] = func(x, y);
    }
  }

  const Z hpitch = n2 * NVAR * sizeof(R); // no ghost zone in the output
  const Z dpitch = s         * sizeof(R);
  cudaMemcpy2D(u, dpitch, host, hpitch, hpitch, n1, cudaMemcpyHostToDevice);
}
