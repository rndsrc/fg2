/* Copyright (C) 2012 Simon Candelaresi
   Copyright (C) 2012 NORDITA

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
#include <cuda_runtime.h>
#include "fg2.h"
#include <hdf5.h>
#include <hdf5_hl.h>

#define N_PARAMS 6

void dump(Z i, R T, const char *ext)
{
  char name[64];
  snprintf(name, sizeof(name), "%04d.%s", i, ext);
  hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  using namespace global;

  hsize_t dims[3]={n1,n2,NVAR};

  // make the table with the parameters
 typedef struct Parameters {
    int    n1;
    int    n2;
    int    n_variables;
    int    size_elements;
    int    i;  // time step
    R      T;  // physical time
  } Parameters;

  size_t dst_size =  sizeof( Parameters );
  size_t dst_offset[N_PARAMS] = { HOFFSET( Parameters, n1 ),
                           HOFFSET( Parameters, n2 ),
                           HOFFSET( Parameters, n_variables ),
                           HOFFSET( Parameters, size_elements ),
                           HOFFSET( Parameters, i ),
                           HOFFSET( Parameters, T )};

  const char *field_names[N_PARAMS]  =
  { "n1","n2", "n_variables", "data size", "time step", "time" };
  Parameters size[1] = {{n1, n2, sizeof(S) / sizeof(R), sizeof(R), i, T}};
  hid_t  field_type[N_PARAMS];
  field_type[0] = H5T_NATIVE_INT;
  field_type[1] = H5T_NATIVE_INT;
  field_type[2] = H5T_NATIVE_INT;
  field_type[3] = H5T_NATIVE_INT;
  field_type[4] = H5T_NATIVE_INT;
  field_type[5] = H5T_NATIVE_FLOAT;
  H5TBmake_table("parameters", file, "table", N_PARAMS, 1,
                         dst_size, field_names, dst_offset, field_type,
                         10, NULL, 0, size);

  const Z hpitch = n2 * sizeof(S); // no ghost zone in the output
  const Z dpitch = s  * sizeof(R);
  cudaMemcpy2D(host, hpitch, u, dpitch, hpitch, n1, cudaMemcpyDeviceToHost);

  // prepare the hdf5 file
  H5LTmake_dataset(file,"/dset",3,dims,H5T_NATIVE_FLOAT,host);

  H5Fclose(file);
}

int load(const char *file_name, R *T) {
  using namespace global;

  typedef struct Parameters {
    int    n1;
    int    n2;
    int    n_variables;
    int    size_elements;
    int    i;  // time step
    R      T;  // physical time
  } Parameters;
  size_t dst_size =  sizeof( Parameters );
  Parameters  dst_buf[1];
  size_t dst_offset[N_PARAMS] = { HOFFSET( Parameters, n1 ),
                           HOFFSET( Parameters, n2 ),
                           HOFFSET( Parameters, n_variables ),
                           HOFFSET( Parameters, size_elements ),
                           HOFFSET( Parameters, i ),
                           HOFFSET( Parameters, T )};
  size_t dst_sizes[N_PARAMS] = { sizeof( dst_buf[0].n1),
                          sizeof( dst_buf[0].n2),
                          sizeof( dst_buf[0].n_variables),
                          sizeof( dst_buf[0].size_elements),
                          sizeof( dst_buf[0].i),
                          sizeof( dst_buf[0].T)};

  hid_t file = H5Fopen (file_name, H5F_ACC_RDONLY, H5P_DEFAULT);

  H5TBread_table(file, "table", dst_size, dst_offset, dst_sizes, dst_buf );

  H5LTread_dataset (file, "dset", H5T_NATIVE_FLOAT, host);

  const Z hpitch = n2 * sizeof(S); // no ghost zone in the output
  const Z dpitch = s  * sizeof(R);
  cudaMemcpy2D(u, dpitch, host, hpitch, hpitch, n1, cudaMemcpyHostToDevice);

  H5Fclose(file);

/*  print("n1 = %i; n2 = %i; n_variables = %i; size_elements = %i; t_idx = %i; T = %f\n",
        dst_buf[0].n1, dst_buf[0].n2, dst_buf[0].n_variables, dst_buf[0].size_elements,
        dst_buf[0].i, dst_buf[0].T);*/

  // copy the parameters to code vartiables


   // TODO: if error return -1
  *T = dst_buf[0].T;
  return dst_buf[0].i;
}
