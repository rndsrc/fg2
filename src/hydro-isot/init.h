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
  const R a  = fi / (nus * ki * ki);

  const R d  = 1.0 + 2.0e-6 * ((double)rand() / RAND_MAX - 0.5);
  const R u2 = a * cos(ki * x);

  return (S){log(d), 0.0, u2};
}

static S zeros(R x, R y)
{
  return (S){0.0, 0.0, 0.0};
}

static S (*pick(const char *name))(R, R)
{
  cudaMemcpyFromSymbol(&nus, "para_nus", sizeof(R));
  cudaMemcpyFromSymbol(&fi,  "para_f",   sizeof(R));
  cudaMemcpyFromSymbol(&ki,  "para_n",   sizeof(R)); ki *= 6.2831853071795865;

  if(!strcmp(name, "Kf")) return Kf;

  return zeros; // default
}
