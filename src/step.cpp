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

int step(const E t, const E dt) // 3rd-order ow-storage Runge-Kutta method
{
  const E alpha[] = {0.0, 1.0/3.0, 3.0/4.0};
  const E beta [] = {0.0, -5.0/9.0, -153.0/128.0};
  const E gamma[] = {1.0/3.0, 15.0/16.0, 8.0/15.0};

  for(Z i = 0; i < 3; ++i) {
    using namespace global;
    bcond(u, 0, 0); // apply periodic boundary condition
    kick (v, u, t + dt * alpha[i], beta[i]); // v <- F(u, T) + beta[i] * v
    drift(u, v,     dt * gamma[i]         ); // u <- u + dt * gamma[i] * v
  }
  cudaThreadSynchronize();

  return (int)cudaGetLastError(); // cudaError_t is enum; cudaSuccess == 0
}
