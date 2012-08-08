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

__device__ __constant__ R para_gamma = 5.0 / 3.0; // ratio of specific heats
__device__ __constant__ R para_nus   = 2.0e-4;    // shear viscosity
__device__ __constant__ R para_nub   = 0.0;       // bulk  viscosity
__device__ __constant__ R para_kappa = 5.0e-4;    // thermal conductivity

__device__ __constant__ R para_dd    = 0.0;       // simple density diffusion

static __device__ S eqns(const S *u, const Z i, const Z j, const Z s)
{
  const S d1 = {D1(lnd), D1(u1), D1(u2), D1(lne)}; // 36 FLOP
  const S d2 = {D2(lnd), D2(u1), D2(u2), D2(lne)}; // 36 FLOP

  S dt = {0.0, 0.0, 0.0, 0.0};

  // Advection: 16 FLOP
  {
    const R u1 = u->u1;
    const R u2 = u->u2;

    dt.lnd -= u1 * d1.lnd + u2 * d2.lnd;
    dt.u1  -= u1 * d1.u1  + u2 * d2.u1 ;
    dt.u2  -= u1 * d1.u2  + u2 * d2.u2 ;
    dt.lne -= u1 * d1.lne + u2 * d2.lne;
  }

  // Compressible/pressure effects: 13 FLOP
  {
    const R gamma1 = para_gamma - K(1.0);
    const R temp   = gamma1 * exp(u->lne);
    const R div_u  = d1.u1 + d2.u2;

    dt.lnd -= div_u;
    dt.u1  -= temp * (d1.lnd + d1.lne);
    dt.u2  -= temp * (d2.lnd + d2.lne);
    dt.lne -= div_u * gamma1;


  // Non-ideal effects (depend on density and temperature): 30 FLOP

    const R s11 =  d1.u1 - div_u  / K(3.0);
    const R s12 = (d1.u2 + d2.u1) / K(2.0);
    const R s22 =  d2.u2 - div_u  / K(3.0);
    const R two_nus = K(2.0) * para_nus;

    dt.u1  += two_nus * (s11 * d1.lnd + s12 * d2.lnd);
    dt.u2  += two_nus * (s12 * d1.lnd + s22 * d2.lnd);
    dt.lne += (gamma1 / temp) *
      (two_nus * (s11 * s11 + K(2.0) * s12 * s12 + s22 * s22) +
       para_nub * div_u * div_u);
  }

  // Non-ideal effects (only depend on velocity): 94 FLOP
  {
    const R d11_u1 = D11(u1), d11_u2 = D11(u2);
    const R d12_u1 = D12(u1), d12_u2 = D12(u2);
    const R d22_u1 = D22(u1), d22_u2 = D22(u2);
    const R mixed  = para_nus / K(3.0) + para_nub;

    dt.u1 += para_nus * (d11_u1 + d22_u1) + mixed * (d11_u1 + d12_u2);
    dt.u2 += para_nus * (d11_u2 + d22_u2) + mixed * (d12_u1 + d22_u2);
  }

  // Density diffusion and thermal conductivity: 64 FLOP
  {
    const R d_lnd_2 = d1.lnd * d1.lnd + d2.lnd * d2.lnd;
    const R d_lne_2 = d1.lne * d1.lne + d2.lne * d2.lne;
    const R ed      = (para_gamma - K(1.0)) * para_kappa;

    dt.lnd += para_dd * (D11(lnd) + D22(lnd) + d_lnd_2);
    dt.lne +=      ed * (D11(lne) + D22(lne) + d_lne_2);
  }

  return dt;
}
