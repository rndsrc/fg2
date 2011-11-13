/* Copyright (C) 2011 Chi-kwan Chan
   Copyright (C) 2011 NORDITA

   This file is part of fd2.

   Fd2 is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Fd2 is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with fd2.  If not, see <http://www.gnu.org/licenses/>. */

#define GET(x, i, j) u[(i) * s + (j)].x

#if ORDER == 6
static __device__ const R coef[] = {45.0/60.0, -9.0/60.0, 1.0/60.0};

#define D1(x) (( coef[0] * (GET(x, 1, 0) - GET(x,-1, 0)) \
               + coef[1] * (GET(x, 2, 0) - GET(x,-2, 0)) \
               + coef[2] * (GET(x, 3, 0) - GET(x,-3, 0)) ) / d1) // 9 FLOP

#define D2(x) (( coef[0] * (GET(x, 0, 1) - GET(x, 0,-1)) \
               + coef[1] * (GET(x, 0, 2) - GET(x, 0,-2)) \
               + coef[2] * (GET(x, 0, 3) - GET(x, 0,-3)) ) / d2) // 9 FLOP
#else
#error !!! Requested derivatives are not implemented !!!
#endif
