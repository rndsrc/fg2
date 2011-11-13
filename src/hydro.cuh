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
  const R ld = u->ld, d1_ld = D1(ld), d2_ld = D2(ld);
  const R u1 = u->u1, d1_u1 = D1(u1), d2_u1 = D2(u1);
  const R u2 = u->u2, d1_u2 = D1(u2), d2_u2 = D2(u2);

  return (S){-(u1 * d1_ld + u2 * d2_ld + d1_u1 + d2_u2),
             -(u1 * d1_u1 + u2 * d2_u1 + d1_ld),
             -(u1 * d1_u2 + u2 * d2_u2 + d2_ld)};
}
