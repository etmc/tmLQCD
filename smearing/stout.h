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

struct stout_control
{
  double                rho;
  unsigned int          iterations;
  unsigned int          current_iteration;
  int                   calculate_force_terms;
  gauge_field_array_t   U;
  gauge_field_array_t   Q;
  complex_field_array_t f0;
  complex_field_array_t f1;
  complex_field_array_t f2;
  real_field_array_t    u;
  real_field_array_t    v;
};

stout_control *construct_stout_control(double rho, unsigned int iterations, int calculate_force_terms);
void free_stout_control(stout_control *control);
gauge_field_t get_stout_smeared_field(stout_control *control);

void stout_smear(stout_control *control, gauge_field_t m_field_in);

void stout_links_with_force_terms(gauge_field_t buff_out, gauge_field_t Q, complex_field_t *f0, complex_field_t *f1, complex_field_t *f2,
                                  real_field_t *u, real_field_t *v, double const coeff, gauge_field_t staples, gauge_field_t buff_in)
