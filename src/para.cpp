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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cuda_runtime.h>
#include "fg2.h"

const char *para(const char *in)
{
  static char sym[256] = "para_";
  sym[5] = '\0'; // necessary because sym is static

  while(*in == '-') ++in; // skip leading dashes
  const Z n = strchr(in, '=') - in;
  strncat(sym, in, n);

  const R val = atof(in + n + 1);
  if(cudaMemcpyToSymbol(sym, &val, sizeof(R)) == cudaSuccess) {
    sprintf(sym + 5 + n, " = %g", val);
    return sym + 5;
  } else
    return NULL;
}
