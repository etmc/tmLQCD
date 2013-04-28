#pragma once

#include <su3.h>

#include <buffers/gauge.h>
#include <buffers/adjoint.h>

/* We need a number of indices to do the bookkeeping.
   This will always amount to at most 12 fields, but
   we define some aliases so that we don't have to do
   the mental mapping all the time. */

void generic_staples(su3 *out, unsigned int x, unsigned int mu, gauge_field_t in);
void generic_staples_3d(su3 *out, unsigned int x, unsigned int mu, gauge_field_t in);
void project_traceless_antiherm(su3 *in);
void reunitarize(su3 *in);

static inline void cayley_hamilton_exponent(su3 *expA, su3 const *A);
void cayley_hamilton_exponent_with_force_terms(su3 *expA, su3 *B1, su3 *B2, _Complex double *f1, _Complex double *f2, su3 const *A);

void unfold_field(gauge_field_t *target, gauge_field_t const base);
void fold_field(gauge_field_t *target, gauge_field_t const base);
void rnd_gauge_trafo(gauge_field_t * target, gauge_field_t const src);

void print_su3(su3 const *in);
void print_su3adj(su3adj const *in);
void print_config_to_screen(gauge_field_t in);

void calculate_forces_numerically(su3adj *result, int const x, int const mu, int * mnllist, const int no);

#include "utils.inline"
