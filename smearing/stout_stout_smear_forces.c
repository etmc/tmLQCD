#include "stout.ih"

#include "stout_construct_Sigma.static"
#include "stout_construct_Z.static"
#include "stout_add_terms_to_forces.static"

void stout_smear_forces(stout_control *control, adjoint_field_t in)
{
  if (!control->calculate_force_terms)
    fatal_error("Smearing control not set up for calculating forces.", "stout_smear_forces");
  
  gauge_field_t smeared_force = get_gauge_field();
  gauge_field_t Sigma_U = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  adjoint_to_gauge(&smeared_force, in);
      
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
#pragma omp parallel
  for (int iter = control->iterations - 1; iter >= 0; --iter)
  {
    /* Since the U's are fields after full smearing steps, U[iter] == U, U[iter + 1] == V */
    construct_Sigma_V(control->trace[iter], control->U[iter + 1], smeared_force);
    construct_Z_V(control->trace[iter], control->U[iter]);
    add_terms_to_forces_V(smeared_force, control->trace[iter], control->rho, control->U[iter + 1], control->U[iter]);

    construct_Sigma_U(Sigma_U, control->rho, control->U[iter + 1], control->trace[iter]);
    add_terms_to_forces_U(smeared_force, Sigma_U, control->U[iter]);
    /* Barrier unnecessary from implicit barrier of single section in add_terms_to_forces */
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(&control->force_result, smeared_force);
  generic_exchange(&control->force_result, sizeof(su3adj_tuple));

  return_gauge_field(&smeared_force);
  return_gauge_field(&Sigma_U);
}
