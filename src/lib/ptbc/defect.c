#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <ptbc.h>
#include "global.h"


// check if a link lie in defect region
// 1. last slice of along 2. 3d cube with 0 origin
bool is_defect(PTBCDefect *def, int const ix, int const mu) {
  if (ix >= VOLUME) return false;
  
  int const global_dim[4] = {T*g_nproc_t, LX*g_nproc_t, LY*g_nproc_x, LZ*g_nproc_z};

  // check start point is at second last slice and link stays on the slice
  if (g_coord[ix][def->along] == global_dim[def->along] - 2 && mu != def->along){
    // xf = end point of link 
    int xf[4] = {g_coord[ix][0], g_coord[ix][1], g_coord[ix][2], g_coord[ix][3]};
    xf[mu]  = xf[mu] + 1;
    int dim3[3];
    int count = 0;
    for (int d=0; d<4; d++) {
      if (d != def->along) {
        dim3[count] = d;
        count++;
      }
    }

    for (int d=0; d<3; d++) {
      if (g_coord[ix][dim3[d]] >= def->Ld[d] || xf[dim3[d]] >= def->Ld[d]) {
        return false;
      }
    }
    if (g_proc_id==0) printf("defect coordinate %d %d %d %d mu = %d\n", g_coord[ix][0], g_coord[ix][1], g_coord[ix][2], g_coord[ix][3], mu);
    return true;
  }
  else {
    return false;
  }
}

// get multiplying factor of parallel tempering locally (assumne no overlapping defects!)
double get_ptbc_coeff(int const ix, int const mu) {
  int const inst = app()->ptbc.instance_id;
  const PTBCInstance *instance = &(app()->ptbc.instances[inst]);
  
  // loop over defects
  for (int i=0; i<instance->n_coeffs; i++) {
    PTBCDefect *def = instance->defects[i];
    // apply coeff if within defect
    if (is_defect(def, ix, mu)) {
      return instance->coefficients[i];
    }
  }

  return 1.;
}
