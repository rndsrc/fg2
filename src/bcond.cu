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

static void periodic(R *x, R *y, Z off, Z h, Z n)
{
  using namespace global;
  const Z pitch = s * sizeof(R);
  const Z width = n * sizeof(R);
  cudaMemcpy2D(x, pitch, x + off, pitch, width, h, cudaMemcpyDeviceToDevice);
  cudaMemcpy2D(y + off, pitch, y, pitch, width, h, cudaMemcpyDeviceToDevice);
}

void bcond(R *x)
{
  using namespace global;
  x -= HALF * (s + NVAR);
  periodic(x, x + HALF * s,    n1 * s,    HALF, (n2 + ORDER) * NVAR);
  periodic(x, x + HALF * NVAR, n2 * NVAR, (n1 + ORDER), HALF * NVAR);
}
