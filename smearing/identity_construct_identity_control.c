#include "identity.ih"

identity_control *construct_identity_control(int calculate_force_terms)
{
  identity_control *control = (identity_control*)malloc(sizeof(identity_control));

  control->smearing_performed = 0;
  
  control->result = NULL;
  control->force_result = NULL;
  
  return control;
}
