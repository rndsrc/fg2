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
#include <cuda_runtime.h>
#include "fg2.h"

void dump(Z i, const char *ext)
{
  char name[64];
  snprintf(name, sizeof(name), "%04d.%s", i, ext);
  FILE *file = fopen(name, "wb");

  using namespace global;

  const Z size[] = {n1, n2, sizeof(S) / sizeof(R), sizeof(R)};
  fwrite(size, sizeof(Z), 4, file);

  const Z hpitch = n2 * sizeof(S); // no ghost zone in the output
  const Z dpitch = s  * sizeof(R);
  cudaMemcpy2D(host, hpitch, u, dpitch, hpitch, n1, cudaMemcpyDeviceToHost);
  fwrite(host, sizeof(S), n1 * n2, file);

  fclose(file);
}
