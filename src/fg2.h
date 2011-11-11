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

typedef int Z;
#if defined(DOUBLE) || defined(OUBLE) // so -DOUBLE works
typedef double R;
#define K(x) (x)
#else
typedef float R;
#define K(x) (x##f)
#endif

void print(const char *, ...);
void error(const char *, ...);
void usage(const char *);

int solve(R, R, Z, Z);

#endif
