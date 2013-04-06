#include "hex.ih"

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const rho_0, double const rho_1, double const rho_2)
{
  hex_control *control = (hex_control*)malloc(sizeof(hex_control));
  control->rho[0] = rho_0;
  control->rho[1] = rho_1;
  control->rho[2] = rho_2;
  control->iterations = iterations;
  control->smearing_performed = 0;
  
  /* We can keep this quite simple if we don't need any forces anyway. */
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1] = get_gauge_field();
   
    control->V_stage_one       = (su3_outer **)malloc(sizeof(su3_outer *));
    control->V_stage_one[0]    = (su3_outer *)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_one   = (stout_notes_outer **)NULL;
    
    control->V_stage_two       = (su3_outer **)malloc(sizeof(su3_outer *));
    control->V_stage_two[0]    = (su3_outer *)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_two   = (stout_notes_outer **)NULL;

    control->trace_stage_three = (stout_notes_tuple **)NULL;
    
    return control;
  }

  control->U = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  
  control->V_stage_one = malloc(iterations * sizeof(su3_outer*)); 
  control->trace_stage_one = (stout_notes_outer **)malloc(iterations * sizeof(stouts_notes_outer*));
    
  control->V_stage_two = malloc(iterations * sizeof(su3_outer*)); 
  control->trace_stage_two = (stout_notes_outer **)malloc(iterations * sizeof(stouts_notes_outer*));

  control->trace_stage_three = (stout_notes_tuple **)malloc(iterations * sizeof(stouts_notes_outer*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {   
    control->V_stage_one[iter]       = (su3_outer*)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_one[iter]   = (stout_notes_outer *)aalloc(sizeof(stout_notes_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->V_stage_two[iter]       = (su3_outer*)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_two[iter]   = (stout_notes_outer *)aalloc(sizeof(stout_notes_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->trace_stage_three[iter] = (stout_notes_tuple *)aalloc(sizeof(stout_notes_tuple) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->U[iter + 1]             = get_gauge_field();
  }
  
  control->result = control->U[iterations]; /* A shallow copy, just putting the reference in place. */
  control->force_result = get_adjoint_field();
  
  return control;
}
