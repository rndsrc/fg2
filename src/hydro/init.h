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

static S implos(R x, R y)
{
  S s = {0.0, 0.0, 0.0, 0.0};

  x -= 0.5 * global::l1;
  y -= 0.5 * global::l2;

  if(fabs(x) + fabs(y) < 0.25) {
    const R d = 0.125;
    const R p = 0.14;
    s.lnd = log(d);
    s.lne = log(p / d / (poly_gamma - 1.0));
  }

  return s;
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

  if(!strcmp(name, "blast" )) return blast;  // Gaussian blast wave
  if(!strcmp(name, "bow"   )) return bow;    // Bow shock
  if(!strcmp(name, "implos")) return implos; // Bow shock
  if(!strcmp(name, "KH"    )) return KH;     // Kelvin-Helmholtz instability
  if(!strcmp(name, "Sod"   )) return Sod;    // Sod shock tube

  return zeros; // default
}
