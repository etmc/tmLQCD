#include "hex.ih"

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const rho_1, double const rho_2, double const rho_3)
{
  hex_control *control = (hex_control*)malloc(sizeof(hex_control));
  control->rho[0] = rho_1;
  control->rho[1] = rho_2;
  control->rho[2] = rho_3;
  control->iterations = iterations;
  control->smearing_performed = 0;
  
  /* We can keep this quite simple if we don't need any forces anyway. */
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1] = get_gauge_field();
    
    control->U_outer = (su3_outer**)malloc(2 * sizeof(su3_outer*));
    control->U_outer[0] = (su3_outer*)aalloc(VOLUMEPLUSRAND * sizeof(su3_outer));
    control->U_outer[1] = (su3_outer*)aalloc(VOLUMEPLUSRAND * sizeof(su3_outer));
    
    return control;
  }

  control->U = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  control->trace = (hex_notes_tuple**)malloc(iterations * sizeof(hex_notes*)); 
  control->U_outer = (su3_outer**)malloc(2 * iterations * sizeof(su3_outer*));
  control->trace_outer = (hex_notes_tuple**)malloc(iterations * sizeof(hex_notes_outer*));  
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {
    control->U[iter + 1] = get_gauge_field();
    control->trace[iter] = aalloc(VOLUME * sizeof(hex_notes));
    control->U_outer[2 * iter] = aalloc(VOLUMEPLUSRAND * sizeof(hex_notes_outer));
    control->U_outer[2 * iter + 1] = aalloc(VOLUMEPLUSRAND * sizeof(hex_notes_outer));
    control->trace_outer[iter] = aalloc(VOLUME * sizeof(hex_notes_outer));
    control->trace_outer[iter + 1] = aalloc(VOLUME * sizeof(hex_notes_outer));
  }
  
  control->result = control->U[iterations]; /* A shallow copy, just putting the reference in place. */
  control->force_result = get_adjoint_field();
  
  return control;
}
