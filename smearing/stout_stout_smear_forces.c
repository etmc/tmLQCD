#include "stout.ih"

#include "stout_add_stout_terms_to_forces.static"
#include "stout_construct_intermediates.static"

void stout_smear_forces(stout_control *control, adjoint_field_t in)
{
  /* Check sanity of the call */
  if (!control->calculate_force_terms)
    fatal_error("Stout control structure not setup for calculating force terms.", "stout_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Stout smearing not yet performed.", "stout_smear_forces");
  
  gauge_field_t smeared_force = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  adjoint_to_gauge(smeared_force, in);
  
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
  for (unsigned int iter = control->iterations - 1; iter > 0; --iter)
  {
    construct_intermediates(control->trace[iter], control->U[iter]);
    add_stout_terms_to_forces(control->trace[iter]->Sigma, control->rho, control->trace[iter]->Z, control->trace[iter], control->U[iter + 1] /* = V */, control->U[iter], smeared_force);
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(control->force_result, sigma);
}
