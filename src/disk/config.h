/* Copyright (C) 2012 Chi-kwan Chan
   Copyright (C) 2012 NORDITA

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

static void config(void)
{
  using namespace global;

  // Simulate the full pi wedge; make grid cells more-or-less square
  l2 = M_PI;
  const R hd2 = 0.5 * l2 / n2;
  const R hd1 = log(hd2 + sqrt(hd2 * hd2 + 1.0));
  l1 = 2.0 * hd1 * n1;

  // Neumann and reflective boundary conditions
  p1 = 0;
  p2 = 0;

  // Compute floating point operation and bandwidth per step
  const Z m1 = n1 + ORDER;
  const Z m2 = n2 + ORDER;
  flops = 3 * ((n1 * n2) * (471 + NVAR * 2.0)); // assume FMA
  bps   = 3 * ((m1 * m2) * 1.0 +
               (n1 * n2) * 5.0 +
               (m1 + m2) * 2.0 * ORDER) * NVAR * sizeof(R) * 8;

  // Set device constant for kernels
  const R Delta[] = {l1 / n1, l2 / n2};
  cudaMemcpyToSymbol("Delta1", Delta+0, sizeof(R));
  cudaMemcpyToSymbol("Delta2", Delta+1, sizeof(R));
}
