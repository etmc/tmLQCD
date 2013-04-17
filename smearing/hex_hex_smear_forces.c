#include "hex.ih"

#include "hex_construct_Sigma_0.static"
#include "hex_construct_Sigma_1.static"
#include "hex_construct_Sigma_2.static"
#include "hex_construct_Sigma_3.static"
#include "hex_construct_Z_1.static"
#include "hex_construct_Z_2.static"
#include "hex_construct_Z_3.static"
#include "hex_add_terms_to_forces_stage_0.static"
#include "hex_add_terms_to_forces_stage_1.static"
#include "hex_add_terms_to_forces_stage_2.static"
#include "hex_add_terms_to_forces_stage_3.static"

void hex_smear_forces(hex_control *control, adjoint_field_t in)
{
  /* Check sanity of the call */
  if (!control->calculate_force_terms)
    fatal_error("Hex control structure not setup for calculating force terms.", "hex_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Hex smearing not yet performed.", "hex_smear_forces");
  
  gauge_field_t smeared_force = get_gauge_field();
  gauge_field_t Sigma_0 = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  adjoint_to_gauge(&smeared_force, in);
      
  /* The modifications are done backwards. */
#pragma omp parallel
  for (unsigned int iter = control->iterations; iter > 0; --iter)
  {
    /* Since the U's are fields after full smearing steps, U[iter - 1] == U, U[iter] == V */
    construct_Sigma_3(control->trace_stage_3[iter - 1], control->U[iter], smeared_force);
    construct_Z_3(control->trace_stage_3[iter - 1], control->U[iter - 1]);
    add_terms_to_forces_stage_3(smeared_force, control->trace_stage_3[iter - 1], control->U[iter], control->U[iter - 1]);
    
    construct_Sigma_2(control->trace_stage_2[iter - 1], control->alpha[2], control->V_stage_2[iter - 1], control->trace_stage_3[iter - 1]);
    construct_Z_2(control->trace_stage_2[iter - 1], control->U[iter - 1] );
    add_terms_to_forces_stage_2(smeared_force, control->trace_stage_2[iter - 1], control->V_stage_2[iter - 1], control->U[iter - 1]);
    
    construct_Sigma_1(control->trace_stage_1[iter - 1], control->alpha[1], control->V_stage_1[iter - 1], control->trace_stage_2[iter - 1]);
    construct_Z_1(control->trace_stage_1[iter - 1], control->U[iter - 1]);
    add_terms_to_forces_stage_1(smeared_force, control->trace_stage_1[iter - 1], control->V_stage_1[iter - 1], control->U[iter - 1]);
    
    construct_Sigma_0(Sigma_0, control->alpha[0], control->U[iter - 1], control->trace_stage_1[iter - 1]);
    add_terms_to_forces_stage_0(smeared_force, Sigma_0, control->U[iter - 1]);
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  gauge_to_adjoint(&control->force_result, smeared_force);
  generic_exchange(&control->force_result, sizeof(su3adj_tuple));

  return_gauge_field(&Sigma_0);
  return_gauge_field(&smeared_force);
}
