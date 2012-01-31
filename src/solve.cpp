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

static const char rotor[] = "-/|\\";

int solve(E t, const E T, Z i, const Z n)
{
  cudaEvent_t c0, c1;
  cudaEventCreate(&c0);
  cudaEventCreate(&c1);
  banner(" Start Simulation ", '=', '=');

  for(const E Dt = (T - t) / (n - i); i++ < n; dump(i, "raw")) {
    const E target = T - (n - i) * Dt; Z m = 0;

    print("%4d:%7.2f ->%7.2f             ", i, (double)t, (double)target);
    cudaEventRecord(c0, 0);
    while(t < target) {
      const E dt = std::min(target - t, getdt());

      print("\b\b\b\b\b\b\b\b\b\b\b\b%c dt ~ %5.0e", rotor[++m%4], (double)dt);
      if(int err = step(t, dt)) {
        print(" crashed in %d step%s\n", m, m > 1 ? "s" : "");
        banner(" Abort Simulation ", '=', '=');
        error("CUDA ERROR: %s\n", cudaGetErrorString((cudaError_t)err));
      }

      t = std::min(t + dt, target);
    }
    cudaEventRecord(c1, 0);
    cudaEventSynchronize(c1);

    using namespace global;
    float ns; cudaEventElapsedTime(&ns, c0, c1); ns *= 1.0e6 / m;
    print("\b\b\b\b\b\b\b\b\b\b\b\b\b, %d step%s "
          "(%.3fms/step, %.3fGflops, %.3fGbps)\n",
          m, m > 1 ? "s" : "", 1.0e-6 * ns, flops / ns, bps / ns);
  }

  banner(" Done  Simulation ", '=', '=');
  cudaEventDestroy(c1);
  cudaEventDestroy(c0);
  return 0;
}
