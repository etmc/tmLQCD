#pragma once

#pragma once

#include <buffers/gauge.h>

typedef struct
{
  double alpha[3];
  unsigned int iterations;

  /* Flags */
  int             smearing_performed;

  /* Result -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */

  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
} hyp_control;

hyp_control *construct_hyp_control(double const alpha_1, double const alpha_2, double const alpha_3,  unsigned int iterations);
void free_hyp_control(hyp_control *control);

void hyp_smear(hyp_control *control, gauge_field_t in);
