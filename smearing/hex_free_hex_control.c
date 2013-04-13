#include "hex.ih"

void free_hex_control(hex_control *control)
{
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  
  if (!control->calculate_force_terms)
  {
    return_gauge_field(&control->U[1]);
    free(control->U);
    
    afree(control->V_stage_0[0]);
    afree(control->V_stage_1[0]);
    
    free(control->V_stage_0);
    free(control->V_stage_1);
    
    free(control);
    return;
  }

  for (unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    afree(control->trace_stage_0[iter]);
    afree(control->V_stage_0[iter]);
    afree(control->trace_stage_1[iter]);
    afree(control->V_stage_1[iter]);
    afree(control->trace_stage_2[iter]);
    return_gauge_field(&control->U[iter + 1]);
  }

  free(control->trace_stage_0);
  free(control->V_stage_0);
  free(control->trace_stage_1);
  free(control->V_stage_1);
  free(control->trace_stage_2);
  
  afree(control->split_forces);

  free(control->U);  
  return_adjoint_field(&control->force_result);
  
  free(control);
}
