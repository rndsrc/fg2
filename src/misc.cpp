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

#include <cstdio>
#include <cstring>
#include "fg2.h"

void banner(const char *title, char l, char r)
{
  int n = strlen(title);
  int h = (80 - 2 - n) / 2, i;
  for(i = 0; i < h; ++i) putchar(l);
  printf("%s%s", title, (n & 1) ? " " : "");
  for(i = 0; i < h; ++i) putchar(r);
  putchar('\n');
  fflush(stdout);
}
