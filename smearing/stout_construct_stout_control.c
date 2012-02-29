#include "stout.ih"

stout_control *construct_stout_control(double rho, unsigned int iterations, int calculate_force_terms)
{
  stout_control *control = (stout_control*)malloc(sizeof(stout_control));
  control->rho = rho;
  control->iterations = iterations;
  control->calculate_force_terms = calculate_force_terms;
  control->smearing_performed = 0;

  /* We can keep this quite simple if we don't need any forces anyway. */
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1] = get_gauge_field();
    control->result = control->U[1];
    return control;
  }

  control->U  = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  control->trace = (stout_notes_tuple**)malloc(iterations * sizeof(stout_notes_tuple*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {
    control->U[iter + 1] = get_gauge_field();
    control->trace[iter] = malloc(VOLUME * sizeof(stout_notes_tuple));
  }
  
  control->result = control->U[iterations]; /* A shallow copy, just putting the reference in place. */
  control->force_result = get_adjoint_field();
  
  return control;
}
