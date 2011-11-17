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
#include <cstring>
#include <cmath>
#include <cuda_runtime.h>
#include "fg2.h"

static R poly_gamma = 5.0 / 3.0;

static S Gaussian(R x, R y)
{
  R d = 0.9 * exp(-0.5 * (x * x + y * y) / 0.01) + 0.1;
  R e = pow(d, poly_gamma) / d;
  return (S){log(d), 0.0, 0.0, log(e)};
}

static S KH(R x, R y)
{
  R d, u;
  if(fabs(y) < 0.25) {
    d =  2.0;
    u =  0.5 + 0.01 * rand() / RAND_MAX - 0.005;
  } else {
    d =  1.0;
    u = -0.5 + 0.01 * rand() / RAND_MAX - 0.005;
  }
  R e = 2.5 / d / (poly_gamma - 1.0);
  return (S){log(d), u, 0.0, log(e)};
}

static S Sod(R x, R y)
{
  R d, P;
  if(fmod(x - 2 * y + K(1.5), K(1.0)) < 0.5) {
    d = 1.0;
    P = 1.0;
  } else {
    d = 0.125;
    P = 0.1;
  }
  R e = P / d / (poly_gamma - 1.0);
  return (S){log(d), 0.0, 0.0, log(e)};
}

void init(const char *name)
{
  S (*func)(R, R) = Gaussian; // default

  if(!strcmp(name, "KH" )) func = KH;  // Kelvin-Helmholtz instability
  if(!strcmp(name, "Sod")) func = Sod; // Sod shock tube

  cudaMemcpyFromSymbol(&poly_gamma, "para_gamma", sizeof(R));

  using namespace global;
  for(Z i = 0; i < n1; ++i) {
    const R x = (i + 0.5) / n1 - 0.5;
    for(Z j = 0; j < n2; ++j) {
      const R y = (j + 0.5) /n2 - 0.5;
      ((S *)host)[i * n2 + j] = func(x, y);
    }
  }

  const Z hpitch = n2 * NVAR * sizeof(R); // no ghost zone in the output
  const Z dpitch = s         * sizeof(R);
  cudaMemcpy2D(u, dpitch, host, hpitch, hpitch, n1, cudaMemcpyHostToDevice);
}
