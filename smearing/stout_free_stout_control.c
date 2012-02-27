#include "stout.ih"

double          rho;
unsigned int    iterations;
unsigned int    current;
int             calculate_force_terms;

void free_stout_control(stout_control *control)
{
  free(control->scratch);
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
    return_gauge_field(&control->Q[iter]);
    return_gauge_field(&control->B1[iter]);
    return_gauge_field(&control->B2[iter]);
    
    free(control->f1[iter]);
    free(control->f2[iter]);
  }
  
  free(control->scratch);
  free(control->U);
  free(control->Q);
  free(control->B1);
  free(control->B2);
  free(control->f1);
  free(control->f2);
  free(control);
  
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  /* The scratch fields should be initialized and cleared by the using functions, not elsewhere. */
}
