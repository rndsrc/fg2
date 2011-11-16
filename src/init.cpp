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

#include <cmath>
#include <cuda_runtime.h>
#include "fg2.h"

static R poly_gamma = 5.0 / 3.0;

static S Gaussian(R x, R y)
{
  x -= 0.5;
  y -= 0.5;

  const R d = 0.9 * exp(-0.5 * (x * x + y * y) / 0.01) + 0.1;
  const R e = pow(d, poly_gamma) / d;

  return (S){log(d), 0.0, 0.0, log(e)};
}

void init(S (*func)(R, R))
{
  cudaMemcpyFromSymbol(&poly_gamma, "para_gamma", sizeof(R));

  if(!func) func = Gaussian;

  using namespace global;

  for(Z i = 0; i < n1; ++i) {
    const R x = (i + 0.5) / n1;
    for(Z j = 0; j < n2; ++j) {
      const R y = (j + 0.5) /n2;
      ((S *)host)[i * n2 + j] = func(x, y);
    }
  }

  const Z hpitch = n2 * NVAR * sizeof(R); // no ghost zone in the output
  const Z dpitch = s         * sizeof(R);
  cudaMemcpy2D(u, dpitch, host, hpitch, hpitch, n1, cudaMemcpyHostToDevice);
}
