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
   
    control->V_stage_0      = (su3_outer **)malloc(sizeof(su3_outer *));
    control->V_stage_0[0]   = (su3_outer *)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_0  = (stout_notes_outer **)NULL;
    
    control->V_stage_1      = (su3_outer **)malloc(sizeof(su3_outer *));
    control->V_stage_1[0]   = (su3_outer *)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_1  = (stout_notes_outer **)NULL;

    control->trace_stage_2  = (stout_notes_tuple **)NULL;
    control->split_forces   = (su3_outer *)NULL:
    
    return control;
  }

  control->U = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  
  control->V_stage_0 = malloc(iterations * sizeof(su3_outer*)); 
  control->trace_stage_0 = (stout_notes_outer **)malloc(iterations * sizeof(stouts_notes_outer*));
    
  control->V_stage_1 = malloc(iterations * sizeof(su3_outer*)); 
  control->trace_stage_1 = (stout_notes_outer **)malloc(iterations * sizeof(stouts_notes_outer*));

  control->trace_stage_2 = (stout_notes_tuple **)malloc(iterations * sizeof(stouts_notes_outer*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {   
    control->V_stage_0[iter]      = (su3_outer*)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_0[iter]  = (stout_notes_outer *)aalloc(sizeof(stout_notes_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->V_stage_1[iter]      = (su3_outer*)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->trace_stage_1[iter]  = (stout_notes_outer *)aalloc(sizeof(stout_notes_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
    
    control->trace_stage_2[iter]  = (stout_notes_tuple *)aalloc(sizeof(stout_notes_tuple) * (VOLUMEPLUSRAND + g_dbw2rand));
    control->U[iter + 1]          = get_gauge_field();
  }
  
  control->split_forces           = (su3_outer*)aalloc(sizeof(su3_outer) * (VOLUMEPLUSRAND + g_dbw2rand));
  
  control->result = control->U[iterations]; /* A shallow copy, just putting the reference in place. */
  control->force_result = get_adjoint_field();
  
  return control;
}
