#include "hex.ih"

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const coeff_1, double const coeff_2, double const coeff_3)
{
  hex_control *control = (hex_control*)malloc(sizeof(hex_control));
  control->coeff[0] = coeff_1;
  control->coeff[1] = coeff_2;
  control->coeff[2] = coeff_3;
  control->iterations = iterations;
  control->calculate_force_terms = calculate_force_terms;
  
  control->result = NULL; /* A shallow copy, just putting the reference in place. */
  
  /* We can keep this quite simple if we don't need any forces anyway. */
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1] = get_gauge_field();
   
    control->V_stage_1      = (su3_two_index **)malloc(sizeof(su3_two_index *));
    control->V_stage_1[0]   = (su3_two_index *)aalloc(sizeof(su3_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_1  = (stout_notes_three_index **)NULL;
    
    control->V_stage_2      = (su3_two_index **)malloc(sizeof(su3_two_index *));
    control->V_stage_2[0]   = (su3_two_index *)aalloc(sizeof(su3_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_2  = (stout_notes_two_index **)NULL;

    control->trace_stage_3  = (stout_notes_tuple **)NULL;
    
    return control;
  }

  control->U = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  
  control->V_stage_1 = malloc(iterations * sizeof(su3_two_index*)); 
  control->trace_stage_1 = (stout_notes_three_index **)malloc(iterations * sizeof(stout_notes_three_index*));
    
  control->V_stage_2 = malloc(iterations * sizeof(su3_two_index*)); 
  control->trace_stage_2 = (stout_notes_two_index **)malloc(iterations * sizeof(stout_notes_two_index*));

  control->trace_stage_3 = (stout_notes_tuple **)malloc(iterations * sizeof(stout_notes_tuple*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {   
    control->V_stage_1[iter]      = (su3_two_index*)aalloc(sizeof(su3_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_1[iter]  = (stout_notes_three_index *)aalloc(sizeof(stout_notes_three_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->V_stage_2[iter]      = (su3_two_index*)aalloc(sizeof(su3_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_2[iter]  = (stout_notes_two_index *)aalloc(sizeof(stout_notes_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->trace_stage_3[iter]  = (stout_notes_tuple *)aalloc(sizeof(stout_notes_tuple) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->U[iter + 1]          = get_gauge_field();
  }
  
  control->force_result = get_adjoint_field();
  
  return control;
}
