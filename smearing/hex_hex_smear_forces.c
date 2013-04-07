#include "hex.ih"

/* The final step of HEX smearing is identical to stouting. */
#include "stout_construct_intermediates.static" 
#include "stout_add_stout_terms_to_forces.static"

#include "hex_construct_intermediates_stage_1.static"
#include "hex_add_hex_terms_to_forces_stage_1.static"

#include "hex_construct_intermediates_stage_0.static"
#include "hex_add_hex_terms_to_forces_stage_0.static"

void hex_smear_forces(hex_control *control, adjoint_field_t in)
{
  /* Check sanity of the call */
  if (!control->calculate_force_terms)
    fatal_error("Hex control structure not setup for calculating force terms.", "hex_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Hex smearing not yet performed.", "hex_smear_forces");
  
  gauge_field_t smeared_force = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  adjoint_to_gauge(&smeared_force, in);
      
  /* The modifications are done backwards. */
#pragma omp parallel
  for (unsigned int iter = control->iterations; iter > 0; --iter)
  {
    stout_construct_intermediates(control->trace_stage_2[iter - 1], control->U[iter] /* = V_3 */, control->U[iter - 1] /* = U */, smeared_force);
    add_hex_terms_to_forces_stage_2(smeared_force, control->rho[2], control->trace_stage_2[iter - 1], control->U[iter] /* = V_3 */, control->U[iter - 1] /* = U */);
    /* What do these forces look like? Shouldn't smeared_force now have 12 components per site, too?*/
    /* Yes, it should. This will be a consequence of adding the staple derivatives! */
    /* So... add_stout_terms_to_forces can't be correct here, can it? It produces a gauge_field_t. */

    construct_intermediates_stage_1(control->trace_stage_1[iter - 1], control->V_stage_1[iter - 1] /* = V_2 */, control->U[iter - 1] /* = U */, smeared_force);
    add_hex_terms_to_forces_stage_1(smeared_force, control->rho, control->trace[iter - 1], control->U[iter] /* = V */, control->U[iter - 1] /* = U */);

    construct_intermediates_stage_0(control->trace_stage_0[iter - 1], control->V_stage_0[iter - 1] /* = V_1 */, control->U[iter - 1] /* = U */, smeared_force);
    add_hex_terms_to_forces_stage_0(smeared_force, control->rho, control->trace[iter - 1], control->U[iter] /* = V */, control->U[iter - 1] /* = U */);
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(&control->force_result, smeared_force);
  generic_exchange(&control->force_result, sizeof(su3adj_tuple));

  return_gauge_field(&smeared_force);
}
