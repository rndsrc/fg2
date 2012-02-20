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

static R M;
static R rS;
static R rin;
static R Gamma;

static S ad_hoc(R lnr, R theta)
{
  const R r   = rin * exp(lnr);
  const R Omg = sin(theta) * sqrt(M / r) / r;

  return (S){0.0, 0.0, 0.0, Omg, 0.0};
}

static S Hawley(R lnr, R theta)
{
  const R r     = rin * exp(lnr);
  const R cyl_r = r * sin(theta);

  // Setup parameters
  const R q0 =  2.0;
  const R r0 = 16.0;
  const R d0 =  1.0;
  const R e0 =  0.01;
  const R bg =  0.0001;

  const R g1 = Gamma - 1.0;
  const R q1 = 2.0 * q0 - 2.0;
  const R lK = sqrt(M * pow(r0, q1 - 1.0));
  const R K  = (g1 * e0) / pow(d0, g1);

  R den = Gamma * e0 - M / r0 + lK * lK / (pow(r0,    q1) * q1)
                      + M / r  - lK * lK / (pow(cyl_r, q1) * q1);
  if(den > 0.0) den =  pow( den * g1 / (Gamma * K), 1.0 / g1);
  else          den = -pow(-den * g1 / (Gamma * K), 1.0 / g1);

  R Omg = lK * pow(cyl_r, -q0);
  if(den < bg) {
    den  = bg;
    Omg *= sin(theta) * sin(theta);
  }
  R eng =  K * pow(den, g1) / g1;

  return (S){log(den), 0.0, 0.0, Omg, log(eng)};
}

static S (*pick(const char *name))(R, R)
{
  cudaMemcpyFromSymbol(&M,     "para_M",     sizeof(R));
  cudaMemcpyFromSymbol(&rS,    "para_rS",    sizeof(R));
  cudaMemcpyFromSymbol(&rin,   "para_rin",   sizeof(R));
  cudaMemcpyFromSymbol(&Gamma, "para_gamma", sizeof(R));

  if(!strcmp(name, "Hawley")) return Hawley; // hydrostatic Hawley (2000) torus

  return ad_hoc; // default
}
