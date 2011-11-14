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

__device__ __constant__ R para_nus   = 5e-4; // shear viscosity
__device__ __constant__ R para_nub   = 0.0;  // bulk  viscosity
__device__ __constant__ R para_gamma = 5.0 / 3.0;

static __device__ S eqns(const S *u, const R d1, const R d2, const Z s)
{
  const R             d1_ld = D1(ld), d2_ld = D2(ld); // 18 FLOP
  const R u1 = u->u1, d1_u1 = D1(u1), d2_u1 = D2(u1); // 18 FLOP
  const R u2 = u->u2, d1_u2 = D1(u2), d2_u2 = D2(u2); // 18 FLOP
  const R             d1_le = D1(le), d2_le = D2(le); // 18 FLOP

  const R s11 = (K(2.0) * d1_u1 - d2_u2) / K(3.0); // 3 FLOP
  const R s12 = (         d1_u2 + d2_u1) / K(2.0); // 2 FLOP
  const R s22 = (K(2.0) * d2_u2 - d1_u1) / K(3.0); // 3 FLOP

  const R gamma_1 = para_gamma - K(1.0);
  const R a = para_nus, b = a / K(3.0) + para_nub, c = a + b; // 3 FLOP

  const R temp = gamma_1 * exp(u->le);
  const R d1_p = temp * (d1_le + d1_ld);
  const R d2_p = temp * (d2_le + d2_ld); // 6 FLOP + FLOP(exp)

  const R diff1 = (K(2.0) * a) * (s11 * d1_ld + s12 * d2_ld)
                +  c * D11(u1) + b * D12(u2) + a * D22(u1); // 49 FLOP
  const R diff2 = (K(2.0) * a) * (s12 * d1_ld + s22 * d2_ld)
                +  a * D11(u2) + b * D21(u1) + c * D22(u2); // 49 FLOP

  const R div_u = d1_u1 + d2_u2;
  const R vis_heating =
    (K(2.0) * para_nus * (s11 * s11 + K(2.0) * s12 * s12 + s22 * s22) +
              para_nub * div_u * div_u) / temp * gamma_1; // 14 FLOP

  return (S){-(u1 * d1_ld + u2 * d2_ld + div_u),
             -(u1 * d1_u1 + u2 * d2_u1 + d1_p - diff1),
             -(u1 * d1_u2 + u2 * d2_u2 + d2_p - diff2),
             -(u1 * d1_le + u2 * d2_le + div_u * gamma_1 - vis_heating)};
                                                                     // 24 FLOP
}
