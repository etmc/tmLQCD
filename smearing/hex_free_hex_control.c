#include "hex.ih"

void free_hex_control(hex_control *control)
{
  if (!control->calculate_force_terms)
  {
    return_gauge_field(&control->U[1]);
    free(control->U);
    afree(control->U_outer[0]);
    afree(control->U_outer[1]);
    free(control);
    return;
  }

  for (unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    return_gauge_field(&control->U[iter + 1]);
    afree(control->trace[iter]);
    afree(control->U_outer[2 * iter]);
    afree(control->U_outer[2 * iter + 1]);
    afree(control->trace_outer[2 * iter]);
    afree(control->trace_outer[2 * iter + 1]);
  }

  return_adjoint_field(&control->force_result);
  
  free(control->U);
  free(control->trace);
  free(control->U_outer);
  free(control->trace_outer);
  free(control);
  
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  /* The scratch fields should be initialized and cleared by the using functions, not elsewhere. */
}
