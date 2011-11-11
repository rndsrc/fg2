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
#include <cuda_runtime.h>
#include "fg2.h"

#define NOVAL (i+1 == argc) || (argv[i+1][0] == '-')
#define BREAK if(NOVAL) break
#define PARA(X) case X: if(NOVAL) goto ignore; // guru can write FORTRAN in C++

int main(int argc, char **argv)
{
  const char *input = "default";

  Z n0 = 10, n1 = 1024, n2 = 1024;
  R t  = 1.;

  int d = 0;

  // If "--help" is an argument, print usage and exit
  for(int i = 1; i < argc; ++i)
    if(!strcmp(argv[i], "--help")) usage(NULL);

  // Home made argument parser
  for(int i = 1; i < argc; ++i) {
    // Arguments do not start with '-' are input files
    if(argv[i][0] != '-') input = argv[i];
    // Arguments start with '-' are options
    else switch(argv[i][1]) {
      PARA('d') d  = atoi(argv[++i]); break;
      PARA('n') n0 = atoi(argv[++i]); BREAK;
           n2 = n1 = atoi(argv[++i]); BREAK;
                n2 = atoi(argv[++i]); break;
      PARA('t') t  = atof(argv[++i]); break;
      default : ignore : usage(argv[i]);
    }
  }
  print("2D finite grid code written in CUDA C\n\n");

  // Pick a device
  int count; cudaGetDeviceCount(&count);
  print("  Device %d/%d : ", d, count);
  if(d < count) {
    if(cudaSuccess == cudaSetDevice(d)) {
      cudaDeviceProp dev; cudaGetDeviceProperties(&dev, d);
      print("\"%s\" with %gMiB of memory\n",
            dev.name, dev.totalGlobalMem / 1048576.0);
    } else
      error("fail to pick device, QUIT\n");
  } else
    error("does not exist, QUIT\n");

  // Setup the grid and global variables
  print("  Resolution : %d x %d", n1, n2);
  if(int sz = setup(n1, n2))
    print(" using %gMiB of memory\n", sz / 1048576.0);
  else
    error(", fail to allocate memory, QUIT\n");

  // TODO: setup initial condition or load starting frame from input
  print("  Initialize : \"%s\"\n", input);

  // Really solve the problem
  print("  Time       : %g with %d frame%s\n", t, n0, n0 > 1 ? "s" : "");
  return solve(0.0, t, 0, n0);
}
