#include "control.ih"

smearing_control_t *construct_smearing_control_from_params(smearing_params_t const *params,  int calculate_force_terms)
{
  switch (params->type)
  {
    case Identity:
      return construct_smearing_control(params->type, calculate_force_terms);
    case APE:
    case APE_3D:
      if (calculate_force_terms)
        fatal_error("APE smearing cannot be used for smearing forces (use Stout instead).", "construct_smearing_control");
    case Stout_3D:
      if (calculate_force_terms)
        fatal_error("Spatial stouting for the HMC algorithm has not yet been implemented.", "construct_smearing_control");      
    case Stout:
      return construct_smearing_control(params->type, calculate_force_terms, params->iterations, params->params[0]);
    case HYP:
    case HYP_3D:
      if (calculate_force_terms)
        fatal_error("HYP smearing cannot be used for smearing forces (use HEX instead).", "construct_smearing_control");
    case HEX_3D:
      if (calculate_force_terms)
        fatal_error("Spatial HEX smearing for the HMC algorithm has not yet been implemented.", "construct_smearing_control");            
    case HEX:
      return construct_smearing_control(params->type, calculate_force_terms, params->iterations, params->params[0], params->params[1], params->params[2]);
    default:
      fatal_error("Requested smearing type not implemented.", "construct_smearing_control_from_params");
  }
  return NULL;
}
