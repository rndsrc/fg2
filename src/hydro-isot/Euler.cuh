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

#include <cstdlib>

__device__ __constant__ R rand1 = 0.0;
__device__ __constant__ R rand2 = 0.0;

static __device__ void first_order(S *u, const R dt, const Z i, const Z j)
{
  if(para_f * para_k > K(0.0)) {
    const R n1 = para_k * cos(rand1);
    const R n2 = para_k * sin(rand1);

    const R phi = K(6.2831853071795865) * ((Z)(n1 + K(0.5)) * i * Delta1 +
                                           (Z)(n2 + K(0.5)) * j * Delta2);
    const R amp = sqrt(dt) * (para_f * cos(phi + rand2)) /
                             (para_k * exp(u->lnd));
    u->u1 += amp * n1;
    u->u2 += amp * n2;
  }
}

#define RANDOMIZE {                                     \
    R r1 = 6.2831853071795865 * rand() / RAND_MAX;      \
    R r2 = 6.2831853071795865 * rand() / RAND_MAX;      \
    cudaMemcpyToSymbol("rand1", &r1, sizeof(R));        \
    cudaMemcpyToSymbol("rand2", &r2, sizeof(R));        \
  }
