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
  R ld;      // ln(density)
  R a, v, l; // specific accretion rate, u_theta.z, specific angular momentum
  R le;      // ln(specific thermal_energy)
};

#ifdef KICK_CU ///////////////////////////////////////////////////////////////

__device__ __constant__ R para_M     = 1.0;       // mass of central black hole
__device__ __constant__ R para_rS    = 2.0;       // Schwarzschild radius
__device__ __constant__ R para_gamma = 5.0 / 3.0; // ratio of specific heats

static __device__ S eqns(const S *u, const Z i, const Z j, const Z s)
{
  S dt = {0.0, 0.0, 0.0, 0.0};

  const R r     = PARA_R0 * exp((i + K(0.5)) * Delta1); // 4 FLOP
  const R r2    = r * r;                                // 1 FLOP
  const R theta =               (j + K(0.5)) * Delta2 ; // 2 FLOP
  const R sin_t = sin(theta);                           // 1 FLOP
  const R cos_t = cos(theta);                           // 1 FLOP

  const R ur     = u->a / r2;          // 1 FLOP
  const R utheta = u->v / sin_t;       // 1 FLOP
  const R uphi   = u->l / (r * sin_t); // 2 FLOP

  const S d1      = {D1(ld), D1(a), D1(v), D1(l), D1(le)}; // 45 FLOP
  const S d2      = {D2(ld), D2(a), D2(v), D2(l), D2(le)}; // 45 FLOP
  const R gamma_1 = para_gamma - K(1.0);                   //  1 FLOP
  const R temp    = gamma_1 * exp(u->le);                  //  2 FLOP
  const R div_u   = (d1.a / r2 + d2.v / sin_t) / r;        //  4 FLOP

  // Advection: 25 FLOP
  {
    dt.ld -= (ur * d1.ld + utheta * d2.ld) / r;
    dt.a  -= (ur * d1.a  + utheta * d2.a ) / r;
    dt.v  -= (ur * d1.v  + utheta * d2.v ) / r;
    dt.l  -= (ur * d1.l  + utheta * d2.l ) / r;
    dt.le -= (ur * d1.le + utheta * d2.le) / r;
  }

  // Compressible/pressure effects: 12 FLOP
  {
    dt.ld -= div_u;
    dt.a  -= temp * (d1.ld + d1.le) * r;
    dt.v  -= temp * (d2.ld + d2.le) / r * sin_t;
    dt.le -= div_u * gamma_1;
  }

  // Pseudo force (coordinate effect): 14 FLOP
  {
    const R up2 = utheta * utheta + uphi * uphi;
    dt.a += (        up2 + K(2.0) * ur * ur    ) * r;
    dt.v += (cos_t * up2 -  sin_t * ur * utheta) / r;
  }

  // External force (gravity): 5 FLOP
  {
    const R r_rS = r - para_rS;
    dt.a -= r2 * para_M / (r_rS * r_rS);
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
  flops = 3 * ((n1 * n2) * (166 + NVAR * 2.0)); // assume FMA
  bps   = 3 * ((m1 * m2) * 1.0 +
               (n1 * n2) * 5.0 +
               (m1 + m2) * 2.0 * ORDER) * NVAR * sizeof(R) * 8;

  // Set device constant for kernels
  const R Delta[] = {l1 / n1, l2 / n2};
  cudaMemcpyToSymbol("Delta1", Delta+0, sizeof(R));
  cudaMemcpyToSymbol("Delta2", Delta+1, sizeof(R));
}

#elif defined(INIT_CPP) //////////////////////////////////////////////////////

static R M;
static R rS;
static R Gamma;

static S ad_hoc(R lnr, R theta)
{
  const R r     = PARA_R0 * exp(lnr);
  const R sin_t = sin(theta);
  const R l     = sqrt(M * r * sin_t * sin_t * sin_t);

  return (S){0.0, 0.0, 0.0, l, 0.0};
}

static S Hawley(R lnr, R theta)
{
  const R r = PARA_R0 * exp(lnr);
  const R p = r * sin(theta);

  // We use something very similar to the steady state torus solution
  // given by Hawley (2000) as our initial condition.  The density is
  // given implicitly by equation (7) in the paper:
  //
  //   Gamma K P / (Gamma - 1) rho = C - Psi - lK^2 / (2q - 2) R^(2q - 2)
  //
  // By comparing the dimensions of different terms on the right hand
  // side, it is clear that R must be dimensionless.  Indeed, Hawley
  // choice the Schwarzschild radius rS = 1.  This automatically gives
  // Psi ~ c^2.
  //
  // Using P = K rho^Gamma, the "pressure acceleration" from the left
  // hand side is
  //
  //   grad(LHS) / rho = Gamma K^2 rho^(Gamma - 2) grad(rho)
  //
  // Comparing this term with grad(P) / rho, it seems the extra
  // polytropic constant "K" is a typo.  We will drop it when we
  // construct our initial condition.
  //
  // The integration constant C controls the size of the torus.  In
  // the Hawley (2000) paper, it is solved by fixing the inner edge of
  // the torus.  Nevertheless, this constant physically controls the
  // temperature, which gives raise to the scale height etc.  To see
  // this, the centrifugal support *almost* cancels gravity at the
  // pressure maximum so we left with
  //
  //   Gamma P / (Gamma - 1) rho = Gamma kB T / (Gamma - 1) mu mH ~ C
  //
  // The contant C determines the temperature at the pressure maximum,
  // and vice versa.  In this file, we will use the temperature, or
  // specific thermal energy, at the pressure maximum, tmp0, (and
  // other parameters) to choose C.

  // Setup parameters
  const R q0 =  2.0;
  const R r0 = 16.0;
  const R d0 =  1.0;
  const R e0 =  0.01;
  const R d1 =  0.01;
  const R e1 =  0.01;

  // Shorthands
  const R g1 = Gamma - 1.0;
  const R q1 = 2.0 * q0 - 2.0;

  // "Specific angular momentum" at pressure maximum: taking the
  // derivative of the right hand side of equation (7) in Hawley
  // (2000), we know the following holds at the pressure maximum
  //
  //   G M / (r - rS)^2 = lK^2 / r^(2 q - 1)
  //
  // Therefore, the following formula is exact and it fixes the unit
  // problem

  const R lK = sqrt(M * pow(r0, q1) * r0) / (r0 - rS);

  // We drop the extra polytropic constant in the left hand side.  We
  // also use specific thermal energy to specify the polytropic and
  // integration constant

  const R K    = (g1 * e0) / pow(d0, g1);
  const R c0   = Gamma * e0 - M / (r0 - rS) + lK * lK / (pow(r0, q1) * q1);
        R prof =         c0 + M / (r  - rS) - lK * lK / (pow(p,  q1) * q1);

  if(prof > 0.0) prof =  pow( prof * g1 / (Gamma * K), 1.0 / g1);
  else           prof = -pow(-prof * g1 / (Gamma * K), 1.0 / g1);

  R d, l, e;
  if(prof > d0 * d1) {
    d = prof;
    l = lK * pow(p, 2.0 - q0);
    e = K * pow(d, g1) / g1;
  } else {
    d = d0 * d1;
    l = 0.0;
    e = K * pow(d, g1) / g1 * e1;
  }

  return (S){log(d), 0.0, 0.0, l, log(e)};
}

static S (*pick(const char *name))(R, R)
{
  cudaMemcpyFromSymbol(&M,     "para_M",     sizeof(R));
  cudaMemcpyFromSymbol(&rS,    "para_rS",    sizeof(R));
  cudaMemcpyFromSymbol(&Gamma, "para_gamma", sizeof(R));

  if(!strcmp(name, "Hawley")) return Hawley; // hydrostatic Hawley (2000) torus

  return ad_hoc; // default
}

#endif ///////////////////////////////////////////////////////////////////////
