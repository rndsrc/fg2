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

static R nus, fi, ki;

static S Kf(R x, R y)
{
  const R r1 = (double)rand() / RAND_MAX;
  const R r2 = (double)rand() / RAND_MAX;
  const R q1 = sqrt(-2.0 * log(r1)) * sin(2.0 * M_PI * r2);
  const R q2 = sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2);

  const R a  = fi / (nus * ki * ki);
  const R u1 = 1.0e-9 * q1;
  const R u2 = 1.0e-9 * q2 + a * cos(ki * x);

  return (S){0.0, u1, u2};
}

static S zeros(R x, R y)
{
  return (S){0.0, 0.0, 0.0};
}

static S (*pick(const char *name))(R, R)
{
  cudaMemcpyFromSymbol(&nus, "para_nus", sizeof(R));
  cudaMemcpyFromSymbol(&fi,  "para_f",   sizeof(R));
  cudaMemcpyFromSymbol(&ki,  "para_k",   sizeof(R));

  if(!strcmp(name, "Kf")) return Kf;

  return zeros; // default
}
