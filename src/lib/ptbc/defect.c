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


// check if a link lie in defect region, cuts at pos + 1/2 and pos + Ld + 1/2
// if either or end of the link is in the defect region, return true 
// include both case when link is in the defect and case link is crossing defect boundary
bool is_defect(PTBCDefect *def, int const ix, int const mu) {
  if (ix >= VOLUME) return false;
  int const global_dim[4] = {T*g_nproc_t, LX*g_nproc_t, LY*g_nproc_x, LZ*g_nproc_z};
  int* coords = g_coord[ix];

  // find three non-mu directions
  int dim3[3];
  int count = 0;
  for (int d=0; d<4; d++) {
    if (d != mu) {
      dim3[count] = d;
      count++;
    }
  }
  // check if mu direction is in the defect or crossing the cut
  if (coords[mu] >= def->pos[mu] && coords[mu] + 1 <= def->pos[mu] + def->Ld[mu] + 1) {
    // find three non-mu directions
    int dim3[3];
    int count = 0;
    for (int d=0; d<4; d++) {
      if (d != mu) {
        dim3[count] = d;
        count++;
      }
    }
    // check if the other three directions are within the defect region
    for (int d=0; d<3; d++) {
      if (coords[dim3[d]] > def->pos[dim3[d]] && coords[dim3[d]] <= def->pos[dim3[d]] + def->Ld[dim3[d]]){
        continue;
      }
      else{
        return false;
      }
    }
    return true;
  }
  else
    return false;  
}


/* bool is_defect(PTBCDefect *def, int const ix, int const mu) {
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
    return true;
  }
  else {
    return false;
  }
} */

// get multiplying factor of parallel tempering locally (assumne no overlapping defects!)
double get_ptbc_coeff(int const ix, int const mu) {
  int const inst = app()->ptbc.instance_id;
  const PTBCInstance *instance = &(app()->ptbc.instances[inst]);

  // if instance is not active, return 1
  if (!instance->active) return 1.;
  
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