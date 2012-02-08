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

struct state {
  R lnd;    // ln(density)
  R u1, u2; // velocity
  R lne;    // ln(specific_thermal_energy)
};

#ifdef KICK_CU ///////////////////////////////////////////////////////////////

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
    const R d_lnd_2 = d1.lnd * d1.lnd + d2.lnd + d2.lnd;
    const R d_lne_2 = d1.lne * d1.lne + d2.lne * d2.lne;
    const R ed      = (para_gamma - K(1.0)) * para_kappa;

    dt.lnd += para_dd * (D11(lnd) + D22(lnd) + d_lnd_2);
    dt.lne +=      ed * (D11(lne) + D22(lne) + d_lne_2);
  }

  return dt;
}

#elif defined(MAIN_CPP) //////////////////////////////////////////////////////

static void config(void)
{
  using namespace global;

  // Compute floating point operation and bandwidth per step
  const Z m1 = n1 + ORDER;
  const Z m2 = n2 + ORDER;
  flops = 3 * ((n1 * n2) * (289 + NVAR * 2.0)); // assume FMA
  bps   = 3 * ((m1 * m2) * 1.0 +
               (n1 * n2) * 5.0 +
               (m1 + m2) * 2.0 * ORDER) * NVAR * sizeof(R) * 8;

  // Set device constant for kernels
  const R Delta[] = {l1 / n1, l2 / n2};
  cudaMemcpyToSymbol("Delta1", Delta+0, sizeof(R));
  cudaMemcpyToSymbol("Delta2", Delta+1, sizeof(R));
}

#elif defined(INIT_CPP) //////////////////////////////////////////////////////

static R poly_gamma;

static R Fermi_Dirac(R x, R d)
{
  return 1.0 / (exp(x / d) + 1.0);
}

static S blast(R x, R y)
{
  x -= 0.5 * global::l1;
  y -= 0.5 * global::l2;

  R e = ((x * x + y * y < 0.01) ? 1.0 : 0.01) / (poly_gamma - 1.0);

  return (S){0.0, 0.0, 0.0, log(e)};
}

static S bow(R x, R y) // need to turn on density diffusion
{
  x -= 0.5 * global::l1;
  y -= 0.5 * global::l2;

  R r = sqrt(x * x + y * y);
  R d = 9.9 * Fermi_Dirac(r - 0.01, 0.001) + 0.1;
  R u = (10.0 / d - 1.0) / 99.0;
  R e = 0.01 / d / (poly_gamma - 1.0);

  return (S){log(d), u, 0.0, log(e)};
}

static S KH(R x, R y)
{
  x -= 0.5 * global::l1;
  y -= 0.5 * global::l2;

  R d, u;
  if(fabs(y) < 0.25) {
    d =  2.0;
    u =  0.5 + 0.01 * rand() / RAND_MAX - 0.005;
  } else {
    d =  1.0;
    u = -0.5 + 0.01 * rand() / RAND_MAX - 0.005;
  }
  R e = 2.5 / d / (poly_gamma - 1.0);

  return (S){log(d), u, 0.0, log(e)};
}

static S Sod(R x, R y)
{
  x -= 0.5 * global::l1;
  y -= 0.5 * global::l2;

  R d, P;
  if(fmod(x - 2 * y + K(1.5), K(1.0)) < 0.5) {
    d = 1.0;
    P = 1.0;
  } else {
    d = 0.125;
    P = 0.1;
  }
  R e = P / d / (poly_gamma - 1.0);

  return (S){log(d), 0.0, 0.0, log(e)};
}

static S zeros(R x, R y)
{
  return (S){0.0, 0.0, 0.0, 0.0};
}

static S (*pick(const char *name))(R, R)
{
  cudaMemcpyFromSymbol(&poly_gamma, "para_gamma", sizeof(R));

  if(!strcmp(name, "blast")) return blast; // Gaussian blast wave
  if(!strcmp(name, "bow"  )) return bow;   // Bow shock
  if(!strcmp(name, "KH"   )) return KH;    // Kelvin-Helmholtz instability
  if(!strcmp(name, "Sod"  )) return Sod;   // Sod shock tube

  return zeros; // default
}

#endif ///////////////////////////////////////////////////////////////////////
