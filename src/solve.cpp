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

#include <algorithm>
#include <cuda_runtime.h>
#include "fg2.h"

int solve(R t, R T, Z i, Z n)
{
  const char rotor[] = "-/|\\";

  cudaEvent_t t0, t1;
  cudaEventCreate(&t0);
  cudaEventCreate(&t1);
  banner(" Start Simulation ", '=', '=');

  for(R dT = (T - t) / n; i++ < n; dump(i, "raw")) {
    Z m = 0;

    print("%4d:%7.2f ->%7.2f             ", i, t, T = dT * i);
    cudaEventRecord(t0, 0);
    while(t < T) {
      const R dt = std::min(T - t, K(1.0) / 1024); // TODO: dynamical dt

      print("\b\b\b\b\b\b\b\b\b\b\b\b%c dt ~ %5.0e", rotor[++m%4], dt);
      step(t, dt);

      t = std::min(t + dt, T);
    }
    cudaEventRecord(t1, 0);
    cudaEventSynchronize(t1);

    using namespace global;
    float ns; cudaEventElapsedTime(&ns, t0, t1); ns *= 1e6 / m;
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b, %d step%s "
          "(%.3fms/step, %.3fGflops, %.3fGbps)\n",
          m, m > 1 ? "s" : "", 1e-6 * ns, flops / ns, bps / ns);
  }

  banner(" Done  Simulation ", '=', '=');
  cudaEventDestroy(t1);
  cudaEventDestroy(t0);
  return 0;
}
