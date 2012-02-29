#include "stout.ih"

void stout_smear(stout_control *control, gauge_field_t in)
{
  control->U[0] = in; /* Shallow copy for convenience */
  
  if (!control->calculate_force_terms)
    stout_smear_plain(control, in);
  else
    stout_smear_with_force_terms(control, in);
}
