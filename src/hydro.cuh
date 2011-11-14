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

__device__ __constant__ R para_gamma = 5.0 / 3.0;
__device__ __constant__ R para_nus   = 5.0e-4;  // shear viscosity
__device__ __constant__ R para_nub   = 0.0;     // bulk  viscosity
__device__ __constant__ R para_kappa = 1.25e-3; // thermal conductivity

static __device__ S eqns(const S *u, const R d1, const R d2, const Z s)
{
  const R gamma_1 = para_gamma - K(1.0), temp = gamma_1 * exp(u->le);
  const R a = para_nus, b = a / K(3.0) + para_nub, c = a + b; // 6 FLOP

  const R             d1_ld = D1(ld), d2_ld = D2(ld); // 18 FLOP
  const R u1 = u->u1, d1_u1 = D1(u1), d2_u1 = D2(u1); // 18 FLOP
  const R u2 = u->u2, d1_u2 = D1(u2), d2_u2 = D2(u2); // 18 FLOP
  const R             d1_le = D1(le), d2_le = D2(le); // 18 FLOP

  const R s11 = (K(2.0) * d1_u1 - d2_u2) / K(3.0);
  const R s12 = (         d1_u2 + d2_u1) / K(2.0);
  const R s22 = (K(2.0) * d2_u2 - d1_u1) / K(3.0); // 8 FLOP

  const R div_u  = d1_u1 + d2_u2; // 1 FLOP
  const R u1_src = (K(2.0) * a) * (s11 * d1_ld + s12 * d2_ld)
                 + c * D11(u1) + b * D12(u2) + a * D22(u1)
                 - temp * (d1_le + d1_ld); // 54 FLOP
  const R u2_src = (K(2.0) * a) * (s12 * d1_ld + s22 * d2_ld)
                 + a * D11(u2) + b * D21(u1) + c * D22(u2)
                 - temp * (d2_le + d2_ld); // 54 FLOP
  const R le_src = gamma_1 *
    (para_kappa * (D11(le) + D22(le) + d1_le * d1_le + d2_le * d2_le) +
     (para_nus * (s11 * s11 + K(2.0) * s12 * s12 + s22 * s22) * K(2.0) +
      para_nub * div_u * div_u) / temp - div_u); // 45 FLOP

  return (S){- (u1 * d1_ld + u2 * d2_ld) - div_u,
             - (u1 * d1_u1 + u2 * d2_u1) + u1_src,
             - (u1 * d1_u2 + u2 * d2_u2) + u2_src,
             - (u1 * d1_le + u2 * d2_le) + le_src}; // 20 FLOP
}
