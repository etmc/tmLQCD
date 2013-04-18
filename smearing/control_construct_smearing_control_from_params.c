#include "control.ih"

smearing_control_t *construct_smearing_control_from_params(smearing_params_t const *params,  int calculate_force_terms)
{
  switch (params->type)
  {
    case Identity:
      return(construct_smearing_control(params->type, calculate_force_terms));
    case APE:
      if (calculate_force_terms)
        fatal_error("APE smearing cannot be used for smearing forces (use Stout instead).", "construct_smearing_control");
    case Stout:
      return(construct_smearing_control(params->type, calculate_force_terms, params->iterations, params->params[0]));
    case HYP:
      if (calculate_force_terms)
        fatal_error("HYP smearing cannot be used for smearing forces (use HEX instead).", "construct_smearing_control");
      return(construct_smearing_control(params->type, calculate_force_terms, params->iterations, params->params[0], params->params[1], params->params[2]));
    case HEX:
      return(construct_smearing_control(params->type, calculate_force_terms, params->iterations, params->params[0], params->params[1], params->params[2]));
    default:
      sprintf(err, "Requested smearing type (%d) implemented.", type);
      fatal_error(err, "construct_smearing_control");
  }
  return NULL;
}
