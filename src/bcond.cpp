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

#include <cuda_runtime.h>
#include "fg2.h"

static void shuffle(R *x, R *y, const Z off, const Z h, const Z n)
{
  using namespace global;
  const Z pitch = s * sizeof(R);
  const Z width = n * sizeof(R);
  cudaMemcpy2D(x, pitch, x + off, pitch, width, h, cudaMemcpyDeviceToDevice);
  cudaMemcpy2D(y + off, pitch, y, pitch, width, h, cudaMemcpyDeviceToDevice);
}

static void periodic1(R *x)
{
  using namespace global;
  x -= HALF * (s + NVAR);
  shuffle(x, x + HALF * s, n1 * s, HALF, (n2 + ORDER) * NVAR);
}

static void periodic2(R *x)
{
  using namespace global;
  x -= HALF * (s + NVAR);
  shuffle(x, x + HALF * NVAR, n2 * NVAR, (n1 + ORDER), HALF * NVAR);
}

static void Neumann1(R *x)
{
  using namespace global;
  const Z width = (n2 + ORDER) * NVAR * sizeof(R);
  x -= HALF * NVAR;
  for(Z i = 0; i < HALF; ++i)
    cudaMemcpy(x - (i + 1) * s,
               x + (i    ) * s, width, cudaMemcpyDeviceToDevice);
  x += n1 * s;
  for(Z i = 0; i < HALF; ++i)
    cudaMemcpy(x + (i    ) * s,
               x - (i + 1) * s, width, cudaMemcpyDeviceToDevice);
}

void bcond(R *x, const int p1)
{
  if(p1) periodic1(x);
  else   Neumann1 (x);

  periodic2(x);
}
