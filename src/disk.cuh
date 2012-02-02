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

#define PARA_R0 K(3.0)

struct state {
  R ld;        // ln(density)
  R ur, ut, l; // u_r, u_theta, specific angular momentum
};

#ifdef KICK_CU ///////////////////////////////////////////////////////////////

static __device__ S eqns(const S *u, const Z i, const Z j, const Z s)
{
  S dt = {0.0, 0.0, 0.0, 0.0};

  const R r     = PARA_R0 * exp((i + K(0.5)) * Delta1); // 4 FLOP
  const R theta =               (j + K(0.5)) * Delta2 ; // 2 FLOP
  const R cot_t = cos(theta) / sin(theta);              // 3 FLOP

  const R ur    = u->ur;
  const R ut    = u->ut;
  const R uphi  = u->l / (r * sin(theta)); // 3 FLOP

  const S d1    = {D1(ld), D1(ur), D1(ut), D1(l)};                // 36 FLOP
  const S d2    = {D2(ld), D2(ur), D2(ut), D2(l)};                // 36 FLOP
  const R div_u = (d1.ur + d2.ut + K(2.0) * ur + ut * cot_t) / r; //  6 FLOP

  // Advection: 20 FLOP
  {
    dt.ld -= (ur * d1.ld + ut * d2.ld) / r;
    dt.ur -= (ur * d1.ur + ut * d2.ur) / r;
    dt.ut -= (ur * d1.ut + ut * d2.ut) / r;
    dt.l  -= (ur * d1.l  + ut * d2.l ) / r;
  }

  // Compressible/pressure effects: 1 FLOP
  {
    dt.ld -= div_u;
  }

  // Pseudo force (coordinate effect): 10 FLOP
  {
    const R uphi2 = uphi * uphi;
    dt.ur += (uphi2         + ut * ut) / r;
    dt.ut += (uphi2 * cot_t - ur * ut) / r;
  }

  return dt;
}

#elif defined(MAIN_CPP) //////////////////////////////////////////////////////

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
  flops = 3 * ((n1 * n2) * (121 + NVAR * 2.0)); // assume FMA
  bps   = 3 * ((m1 * m2) * 1.0 +
               (n1 * n2) * 5.0 +
               (m1 + m2) * 2.0 * ORDER) * NVAR * sizeof(R) * 8;

  // Set device constant for kernels
  const R Delta[] = {l1 / n1, l2 / n2};
  cudaMemcpyToSymbol("Delta1", Delta+0, sizeof(R));
  cudaMemcpyToSymbol("Delta2", Delta+1, sizeof(R));
}

#elif defined(INIT_CPP) //////////////////////////////////////////////////////

static S ad_hoc(R r, R theta)
{
  return (S){0.0, 0.0, 0.0, 0.0};
}

static S (*pick(const char *name))(R, R)
{
  return ad_hoc; // default
}

#endif ///////////////////////////////////////////////////////////////////////
