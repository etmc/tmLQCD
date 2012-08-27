#include "stout.ih"

void free_stout_control(stout_control *control)
{
  if (!control)
    return;
  if (!control->calculate_force_terms)
  {
    return_gauge_field(&control->U[1]);
    free(control->U);
    free(control);
    return;
  }

  for (unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    return_gauge_field(&control->U[iter + 1]);
    free(control->trace[iter]);
  }

  return_adjoint_field(&control->force_result);
  
  free(control->U);
  free(control->trace);
  free(control);
  
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  /* The scratch fields should be initialized and cleared by the using functions, not elsewhere. */
}
