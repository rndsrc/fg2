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
#include "fg2.h"

void usage(const char *bad)
{
  if(bad)
    fprintf(stderr, "hydro: illegal option \"%s\"\n\
Try `hydro --help` for more information.\n", bad);
  else
    print("Usage: hydro [OPTION...] [INPUT_FILE]\n\
2D finite grid code written in CUDA C\n\
\n\
      --help        display this help and exit\n\
  -d                specify device id\n\
  -l                total time, box size\n\
  -n                number of frames and grids\n\
  parameter=val     set parameters in a particular scheme\n\
\n\
Report bugs to <ckch@nordita.org>.\n");

  exit(bad ? -1 : 0);
}
