#pragma once

#include <buffers/adjoint.h>
#include <buffers/gauge.h>
#include <buffers/real.h>
#include <buffers/spinor.h>

#include <smearing/utils.h>

/* Since smearing with the prospect of calculating the modification to any forces requires
 * a large amount of intermediate data storage, we need some concise way of storing all
 * this, as well as a mechanism for passing all the necessary parameters. The current setup
 * intends to provide this through a control structure that gives maximum flexibility with
 * minimal interface clutter.
 */

typedef _Complex double exp_par[4];

typedef struct
{
  double          rho;
  unsigned int    iterations;
  unsigned int    current;
  int             calculate_force_terms;
  int             smearing_performed;
  
  gauge_field_t    result; /* For direct access to the result, shallow copy... */
  gauge_field_t   *scratch;
  adjoint_field_t  force_result;
  
  /* The following are fields that store intermediate results for the force terms */
  gauge_field_t   *U;
  gauge_field_t   *Q;
  gauge_field_t   *B1;
  gauge_field_t   *B2;
  
  exp_par        **f1;
  exp_par        **f2;
} stout_control;

stout_control *construct_stout_control(double rho, unsigned int iterations, int calculate_force_terms);
void free_stout_control(stout_control *control);

void stout_smear(stout_control *control, gauge_field_t in);
void stout_smear_forces(stout_control *control, adjoint_field_t in);
