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

static __device__ S eqns(const S *u, const R d1, const R d2, const Z s)
{
  const R dt_den = - u->u1 * (u[0].den - u[-s].den) / d1
                   - u->u2 * (u[0].den - u[-1].den) / d2;

  return (S){dt_den, K(0.0), K(0.0)};
}
