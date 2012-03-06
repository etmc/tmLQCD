#include "stout.ih"

#include "stout_add_stout_terms_to_forces.static"
#include "stout_construct_lambda.static"

void stout_smear_forces(stout_control *control, adjoint_field_t in)
{
  /* Check sanity of the call */
  if (!control->calculate_force_terms)
    fatal_error("Stout control structure not setup for calculating force terms.", "stout_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Stout smearing not yet performed.", "stout_smear_forces");

  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  gauge_field_t sigma = get_gauge_field();
  adjoint_to_gauge(sigma, in);

  /* The staple derivatives contain non-local factors of lambda, so we'll need to store a whole field of those. */
  gauge_field_t lambda = get_gauge_field();

  
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
  for (unsigned int iter = control->iterations - 1; iter > 0; --iter)
  {
    construct_lambda(lambda, sigma, control->trace[iter], control->U[iter]);
    add_stout_terms_to_forces(sigma, control->rho, lambda, control->trace[iter], control->U[iter]);
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(control->force_result, sigma);

  return_gauge_field(&lambda);
  return_gauge_field(&sigma);
}
