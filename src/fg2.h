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

#ifndef FG2_H
#define FG2_H

#ifndef ORDER
#define ORDER 6
#endif

#if (ORDER) < 0 || (ORDER) % 2
#error !!! ORDER must be positive and even !!!
#else
#define HALF ((ORDER) / 2)
#endif

typedef int Z;
#if defined(DOUBLE) || defined(OUBLE) // so -DOUBLE works
typedef double R;
#define K(x) (x)
#else
typedef float R;
#define K(x) (x##f)
#endif

typedef struct state S;
#define NVAR (sizeof(S) / sizeof(R))

namespace global {
  extern Z n1, n2, s;          // resolution, stride in R
  extern R *u, *v, *host;      // state, swap, and host array
  extern Z g1, g2, b1, b2, sz; // grid and block dim, shared memory in byte
}

void print(const char *, ...);
void error(const char *, ...);
void usage(const char *);

void dump(Z, const char *);

int setup(Z, Z);
int solve(R, R, Z, Z);

struct state {
  R den, u1, u2; // density and velocity
};

#endif
