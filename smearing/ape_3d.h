#pragma once

#include <buffers/gauge.h>

typedef struct
{
  double rho;
  unsigned int iterations;
  
  /* Result -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */

  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
} ape_3d_control;

ape_3d_control *construct_ape_3d_control(unsigned int iterations, double rho);
void free_ape_3d_control(ape_3d_control *control);

void ape_3d_smear(ape_3d_control *control, gauge_field_t in);
