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

__device__ __constant__ R para_gamma = 5.0 / 3.0; // ratio of specific heats
__device__ __constant__ R para_dd    = 0.0;       // density diffusion
__device__ __constant__ R para_nus   = 2.0e-4;    // shear viscosity
__device__ __constant__ R para_nub   = 0.0;       // bulk  viscosity
__device__ __constant__ R para_kappa = 5.0e-4;    // thermal conductivity

static __device__ S eqns(const S *u, const Z s)
{
  S dt = {0.0, 0.0, 0.0, 0.0};

  const S d1      = {D1(ld), D1(u1), D1(u2), D1(le)}; // 36 FLOP
  const S d2      = {D2(ld), D2(u1), D2(u2), D2(le)}; // 36 FLOP
  const R gamma_1 = para_gamma - K(1.0);              // 1 FLOP
  const R temp    = gamma_1 * exp(u->le);             // 1 FLOP
  const R div_u   = d1.u1 + d2.u2;                    // 1 FLOP

  // Advection: 16 FLOP
  {
    const R u1 = u->u1;
    const R u2 = u->u2;

    dt.ld -= u1 * d1.ld + u2 * d2.ld;
    dt.u1 -= u1 * d1.u1 + u2 * d2.u1;
    dt.u2 -= u1 * d1.u2 + u2 * d2.u2;
    dt.le -= u1 * d1.le + u2 * d2.le;
  }

  // Compressible/pressure effects: 9 FLOP
  {
    dt.ld -= div_u;
    dt.u1 -= temp * (d1.ld + d1.le);
    dt.u2 -= temp * (d2.ld + d2.le);
    dt.le -= div_u * gamma_1;
  }

  // Non-ideal effects (only depend on velocity): 94 FLOP
  {
    const R d11_u1 = D11(u1), d11_u2 = D11(u2);
    const R d12_u1 = D12(u1), d12_u2 = D12(u2);
    const R d22_u1 = D22(u1), d22_u2 = D22(u2);
    const R mixed = para_nus / K(3.0) + para_nub;

    dt.u1 += para_nus * (d11_u1 + d22_u1) + mixed * (d11_u1 + d12_u2);
    dt.u2 += para_nus * (d11_u2 + d22_u2) + mixed * (d12_u1 + d22_u2);
  }

  // Non-ideal effects (also depend on density and temperature): 30 FLOP
  {
    const R s11 =  d1.u1 - div_u  / K(3.0);
    const R s12 = (d1.u2 + d2.u1) / K(2.0);
    const R s22 =  d2.u2 - div_u  / K(3.0);
    const R two_nus = K(2.0) * para_nus;

    dt.u1 += two_nus * (s11 * d1.ld + s12 * d2.ld);
    dt.u2 += two_nus * (s12 * d1.ld + s22 * d2.ld);
    dt.le += (gamma_1 / temp) *
      (two_nus * (s11 * s11 + K(2.0) * s12 * s12 + s22 * s22) +
       para_nub * div_u * div_u);
  }

  // Density diffusion and thermal conductivity: 63 FLOP
  {
    const R ed = gamma_1 * para_kappa;
    dt.ld += para_dd * (D11(ld) + D22(ld) + d1.ld * d1.ld + d2.ld + d2.ld);
    dt.le +=      ed * (D11(le) + D22(le) + d1.le * d1.le + d2.le * d2.le);
  }

  return dt;
}
