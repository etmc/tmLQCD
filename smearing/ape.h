#pragma once

#include <buffers/gauge.h>

typedef struct
{
  double coeff;
  unsigned int iterations;
  
  /* Result -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */

  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t   *U;     /* The sequence of iterations gauge fields */
} ape_control;

ape_control *construct_ape_control(unsigned int iterations, double coeff);
void free_ape_control(ape_control *control);

void ape_smear(ape_control *control, gauge_field_t in);
