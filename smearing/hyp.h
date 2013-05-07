#pragma once

#include <buffers/gauge.h>

typedef su3 su3_two_index[12];

typedef struct
{
  double coef[3];
  unsigned int iterations;

  /* Result -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */

  /* Intermediate results, stored to enhance locality of the analysis */
  su3_two_index     **staples; /* Scratch space, available here so it is persistent */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
} hyp_control;

hyp_control *construct_hyp_control(unsigned int iterations, double const coef_0, double const coef_1, double const coef_2);
void free_hyp_control(hyp_control *control);

void hyp_smear(hyp_control *control, gauge_field_t in);
