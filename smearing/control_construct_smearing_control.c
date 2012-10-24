#include "control.ih"

smearing_control *construct_smearing_control(smearing_type type, ...)
{
  /* Declare some local variables to fill from the variadic arguments */
  double params_double[3];
  int params_int[1];
  unsigned int params_uint[1];

  smearing_control *result = (smearing_control*) malloc(sizeof(smearing_control));
  result->type = type;

  va_list smearing_args;
  va_start(smearing_args, type);
  switch (type)
  {
   case Identity:
      result->control = (void*)construct_identity_control();
      break;
   case Stout:
      params_double[0] = va_arg(smearing_args, double);
      params_uint[0] = va_arg(smearing_args, unsigned int);
      params_int[0] = va_arg(smearing_args, int);
      result->control = (void*)construct_stout_control(params_double[0], params_uint[0], params_int[0]);
      break;
  }
  va_end(smearing_args);

  return result;
}

