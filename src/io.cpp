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

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "fg2.h"

static Z frame(const char *h)
{
  char c;
  while((int)(c = *h++))
    if('0' <= c && c <= '9')
      return atoi(h-1); // get the frame number
  return 0;
}

Z load(const char *name)
{
  FILE *file = fopen(name, "rb"); // assume checked with exist()

  using namespace global;

  Z size[4];
  fread(size, sizeof(Z), 4, file);
  if(size[0] != n1 || size[1] != n2 || size[2] != NVAR || size[3] != sizeof(R))
    error("input data is not compatible with the setup");

  const Z hpitch = n2 * NVAR * sizeof(R); // no ghost zone in the output
  const Z dpitch = s         * sizeof(R);
  fread(host, NVAR * sizeof(R), n1 * n2, file);
  cudaMemcpy2D(u, dpitch, host, hpitch, hpitch, n1, cudaMemcpyHostToDevice);

  fclose(file);

  return frame(name);
}

void dump(const Z i, const char *ext)
{
  char name[64];
  snprintf(name, sizeof(name), "%04d.%s", i, ext);
  FILE *file = fopen(name, "wb");

  using namespace global;

  const Z size[] = {n1, n2, NVAR, sizeof(R)};
  fwrite(size, sizeof(Z), 4, file);

  const Z hpitch = n2 * NVAR * sizeof(R); // no ghost zone in the output
  const Z dpitch = s         * sizeof(R);
  cudaMemcpy2D(host, hpitch, u, dpitch, hpitch, n1, cudaMemcpyDeviceToHost);
  fwrite(host, NVAR * sizeof(R), n1 * n2, file);

  fclose(file);
}
