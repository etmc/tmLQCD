#ifndef GUARD_SMEARING_STOUT_3D_H
#define GUARD_SMEARING_STOUT_3D_H

#include <buffers/adjoint.h>
#include <buffers/gauge.h>
#include <buffers/utils.h>

#include <smearing/utils.h>

typedef struct
{
  /* Flags */
  int calculate_force_terms;
  
  /* Parameters */
  double          rho; /* For now, we're going to work with homogeneous smearing coefficients */
  unsigned int    iterations;
    
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
  
  /* Final results -- the first is a shallow copy */
  gauge_field_t    result;
} stout_3d_control;

stout_3d_control *construct_stout_3d_control(int calculate_force_terms, unsigned int iterations, double rho);
void free_stout_3d_control(stout_3d_control *control);

void stout_3d_smear(stout_3d_control *control, gauge_field_t in);

#endif
