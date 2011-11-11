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
#include <cstring>
#include "fg2.h"

#define NOVAL (i+1 == argc) || (argv[i+1][0] == '-')
#define BREAK if(NOVAL) break
#define PARA(X) case X: if(NOVAL) goto ignore; // guru can write FORTRAN in C++

int main(int argc, char **argv)
{
  const char *input = "default";

  Z n0 = 10, n1 = 1024, n2 = 1024;

  // If "--help" is an argument, print usage and exit
  for(int i = 1; i < argc; ++i)
    if(!strcmp(argv[i], "--help")) usage(NULL);

  // Home made argument parser
  for(int i = 1; i < argc; ++i) {
    // Arguments do not start with '-' are input files
    if(argv[i][0] != '-') input = argv[i];
    // Arguments start with '-' are options
    else switch(argv[i][1]) {
      PARA('n') n0 = atoi(argv[++i]); BREAK;
           n2 = n1 = atoi(argv[++i]); BREAK;
                n2 = atoi(argv[++i]); break;
      default : ignore : usage(argv[i]);
    }
  }

  // Print simulation setup
  print("2D finite grid code written in CUDA C\n\n");
  print("  Resolution : %d x %d\n", n1, n2);
  print("  Initialize : \"%s\"\n", input);
  print("  Time       : %d\n", n0);

  return 0;
}
