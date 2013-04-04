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
  adjoint_to_gauge(&smeared_force, in);
      
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
#pragma omp parallel
  for (unsigned int iter = control->iterations; iter > 0; --iter)
  {
    construct_intermediates(control->trace[iter - 1], control->U[iter] /* = V */, control->U[iter - 1] /* = U */, smeared_force);
    add_stout_terms_to_forces(smeared_force, control->rho, control->trace[iter - 1], control->U[iter] /* = V */, control->U[iter - 1] /* = U */);
    /* Barrier unnecessary from implicit barrier of critical section in add_stout_terms_to_forces */
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(&control->force_result, smeared_force);
  generic_exchange(&control->force_result, sizeof(su3adj_tuple));

  return_gauge_field(&smeared_force);
}
