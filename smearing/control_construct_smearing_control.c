#include "control.ih"

smearing_control_t *construct_smearing_control(smearing_type type, int calculate_force_terms, ...)
{
  /* Declare some local variables to fill from the variadic arguments */
  double params_double[3];
  unsigned int params_uint[1];

  smearing_control_t *result = (smearing_control_t*)malloc(sizeof(smearing_control_t));
  result->type = type;
  result->id = -1;

  result->calculate_force_terms = calculate_force_terms; 
  result->smearing_performed = 0;
 
  va_list smearing_args;
  va_start(smearing_args, calculate_force_terms);
  switch (type)
  {
    case Identity:
      result->type_control = (void*)construct_identity_control(calculate_force_terms);
      break;
    case APE:
      if (calculate_force_terms)
        fatal_error("APE smearing cannot be used for smearing forces (use Stout instead).", "construct_smearing_control");
      params_uint[0] = va_arg(smearing_args, unsigned int);
      params_double[0] = va_arg(smearing_args, double);
      result->type_control = (void*)construct_ape_control(params_uint[0], params_double[0]);
      break;
    case HYP:
      if (calculate_force_terms)
        fatal_error("HYP smearing cannot be used for smearing forces (use HEX instead).", "construct_smearing_control");
      params_uint[0] = va_arg(smearing_args, unsigned int);
      params_double[0] = va_arg(smearing_args, double);
      params_double[1] = va_arg(smearing_args, double);
      params_double[2] = va_arg(smearing_args, double);
      result->type_control = (void*)construct_hyp_control(params_uint[0], params_double[0], params_double[1], params_double[2]);
      break;
    case Stout:
      params_uint[0] = va_arg(smearing_args, unsigned int);
      params_double[0] = va_arg(smearing_args, double);
      result->type_control = (void*)construct_stout_control(calculate_force_terms, params_uint[0], params_double[0]);
      break;
    default:
      fatal_error("Requested smearing type not implemented.", "construct_smearing_control");
  }
  va_end(smearing_args);

  result->result = NULL;
  result->force_result = NULL;
  
  return result;
}
