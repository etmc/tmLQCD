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


/* All defined in terms of arrays of tuples -- needed to allow for g_gauge_field as input */

void hyp_staples_exclude_none(gauge_field_t buff_out, gauge_field_array_t buff_in); /* 12 components in, 4 components out */
void hyp_staples_exclude_one (gauge_field_array_t buff_out, gauge_field_array_t buff_in);  /* 12 components in, 12 components out */
void hyp_staples_exclude_two (gauge_field_array_t buff_out, gauge_field_t buff_in);  /*  4 components in, 12 components out */

void APE_project_exclude_one (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
void APE_project_exclude_two (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
