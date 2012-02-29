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

typedef struct
{
  su3    Q;
  su3    B1;
  su3    B2;
  double f1;
  double f2;
} stout_notes_t;

typedef stout_notes_t stout_notes_tuple[4];

typedef struct
{
  /* Parameters */
  double          rho; /* For now, we're going to work with homogeneous smearing coefficients */
  unsigned int    iterations;
  
  /* Flags / trackers */
  unsigned int    current;
  int             calculate_force_terms;
  int             smearing_performed;
  
  /* Results -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */
  adjoint_field_t  force_result;
  
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
  stout_notes_tuple **trace; /* Intermediate results to avoid double calculations */
} stout_control;

stout_control *construct_stout_control(double rho, unsigned int iterations, int calculate_force_terms);
void free_stout_control(stout_control *control);

void stout_smear(stout_control *control, gauge_field_t in);
void stout_smear_forces(stout_control *control, adjoint_field_t in);
