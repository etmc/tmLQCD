#include "stout.ih"

stout_control *construct_stout_control(int calculate_force_terms, unsigned int iterations, double coef)
{
  stout_control *control = (stout_control*)malloc(sizeof(stout_control));
  control->coef = coef;
  control->iterations = iterations;
  
  control->calculate_force_terms = calculate_force_terms;
  control->result = NULL; /* Will be set once the smearing has been performed */

  /* We can keep this quite simple if we don't need any forces anyway. */
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1] = get_gauge_field_annotated("STOUT CONTROL U[1]");
    printf("Field being cleared: id -- %d, bytes -- %llu, note -- %s.\n", ameta_id(control->U[1]), ameta_bytes(control->U[1]), ameta_note(control->U[1]));
    return control;
  }

  control->U  = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  control->trace = (stout_notes_tuple**)malloc(iterations * sizeof(stout_notes_tuple*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {
    control->U[iter + 1] = get_gauge_field();
    control->trace[iter] = aalloc(sizeof(stout_notes_tuple) * (VOLUMEPLUSRAND + g_dbw2rand) + 1);
  }
  control->force_result = get_adjoint_field();
  
  return control;
}
