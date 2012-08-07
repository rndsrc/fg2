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
#include <cstring>
#include <cuda_runtime.h>
#include "fg2.h"

#define ISAB(X) (('A' <= (X) && (X) <= 'Z') || ('a' <= (X) && (X) <= 'z'))
#define NOVAL   (i+1 == argc) || (argv[i+1][0] == '-' && ISAB(argv[i+1][1]))
#define BREAK   if(NOVAL) break
#define PARA(X) case X: if(NOVAL) goto ignore; // guru can write FORTRAN in C++

int main(int argc, char **argv)
{
  const char *input = "default";

  Z i, n0 = 10, n1 = 1024, n2 = 1024, d = 0;
  R t, l0 = 1.0;

  // If "--help" is an argument, print usage and exit
  for(i = 1; i < argc; ++i)
    if(!strcmp(argv[i], "--help")) usage(NULL);

  // Home made argument parser
  for(i = 1; i < argc; ++i) {
    // Check for parameters, we can only set them after picking device
    if(strchr(argv[i], '='));
    // Arguments do not start with '-' are input files
    else if(argv[i][0] != '-') input = argv[i];
    // Arguments start with '-' are options
    else switch(argv[i][1]) {
      PARA('b') global::p1 = atop(argv[++i]); BREAK;
                global::p2 = atop(argv[++i]); break;
      PARA('c') global::c  = atof(argv[++i]); break;
      PARA('d')         d  = atoi(argv[++i]); break;
      PARA('l')         l0 = atof(argv[++i]); BREAK;
   global::l2 = global::l1 = atof(argv[++i]); BREAK;
                global::l2 = atof(argv[++i]); break;
      PARA('n')         n0 = atoi(argv[++i]); BREAK;
                   n2 = n1 = atoi(argv[++i]); BREAK;
                        n2 = atoi(argv[++i]); break;
      default : ignore : usage(argv[i]);
    }
  }

  // TODO: print version and compile options
  print("2D finite grid code written in CUDA C\n\n");

  // Pick a device, obtain global and shared memory size
  double gsz = 0.0, ssz = 0.0;
  cudaGetDeviceCount(&i);
  print("  Device %d/%d   : ", d, i);
  if(d < i) {
    if(cudaSuccess == cudaSetDevice(d)) {
      cudaDeviceProp dev; cudaGetDeviceProperties(&dev, d);
      gsz = dev.totalGlobalMem;
      ssz = dev.sharedMemPerBlock;
      print("\"%s\" with %gMiB global and %gKiB shared memory\n",
            dev.name, gsz / 1048576.0, ssz / 1024.0);
    } else
      error("fail to pick device, QUIT\n");
  } else
    error("does not exist, QUIT\n");

  // Setup the grid and global variables
  print("  Resolution   : %d x %d", n1, n2);
  if(size_t sz = setup(n1, n2))
    print(" using %.3gMiB (%.3g%%) of global memory\n",
          sz / 1048576.0, 100.0 * sz / gsz);
  else
    error(", fail to allocate memory, QUIT\n");

  print("  Grid x block : (%d x %d) x (%d x %d)",
        global::g1, global::g2, global::b1, global::b2);
  print(" using %.3gKiB (%.3g%%) of shared memory\n",
        global::sz / 1024.0, 100.0 * global::sz / ssz);

  // Set parameters
  for(i = 1; i < argc; ++i)
    if(strchr(argv[i], '=')) {
      if(const char *in = para(argv[i]))
        print("  Parameter    : %s\n", in);
      else
        print("  Fail to set  : \"%s\"\n", argv[i]);
    }

  // Setup initial condition or load starting frame from input
  if(exist(input)) {
    print("  Input file   : \"%s\"...", input);
    i = frame(input);
    t = load (input);
    config();
    print("\b\b\b, start at frame %d, time %g\n", i, t);
  } else {
    print("  Initialize   : \"%s\"\n", input);
    config();
    init(input);
    dump(name(i = 0), t = 0.0);
  }

  // Really solve the problem
  print("  Time         : %g with %d frame%s\n", l0, n0, n0 > 1 ? "s" : "");
  return solve(t, t + l0, i, i + n0);
}
