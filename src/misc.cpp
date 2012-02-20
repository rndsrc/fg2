/* Copyright (C) 2011,2012 Chi-kwan Chan
   Copyright (C) 2011,2012 NORDITA

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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "fg2.h"

int atop(const char *s)
{
  return s[0] == 'p' || s[0] == 'P';
}

void banner(const char *title, const char l, const char r)
{
  const Z n = strlen(title);
  const Z h = (80 - 2 - n) / 2;
  for(Z i = 0; i < h; ++i) putchar(l);
  printf("%s%s", title, (n & 1) ? " " : "");
  for(Z i = 0; i < h; ++i) putchar(r);
  putchar('\n');
  fflush(stdout);
}

int exist(const char *name)
{
  FILE *file = fopen(name, "r");
  if(!file) return 0;
  fclose(file);
  return 1;
}

Z frame(const char *name)
{
  char c;
  while((int)(c = *name++))
    if('0' <= c && c <= '9')
      return atoi(name - 1); // get the frame number
  return 0;
}

const char *name(Z i)
{
  static char name[64];
  snprintf(name, sizeof(name), "%04d.%s", i, "raw");
  return name;
}
