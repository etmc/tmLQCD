#pragma once

#include <smearing/stout.h>
#include <smearing/hyp.h>

typedef su3 su3_two_index_3d[6];

typedef struct
{ 
  /* Parameters */
  double          alpha[2];
  unsigned int    iterations;
  
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t  *U; /* The sequence of iterations gauge fields */
  su3_two_index  *V;
  
  /* Final results -- the first is a shallow copy */
  gauge_field_t    result;
} hex_3d_control;

hex_3d_control *construct_hex_3d_control(unsigned int iterations, double const alpha_1, double const alpha_2);
void free_hex_3d_control(hex_3d_control *control);

void hex_3d_smear(hex_3d_control *control, gauge_field_t in);
