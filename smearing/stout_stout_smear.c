#include "stout.ih"

#include "stout_smear_plain.static"
#include "stout_smear_with_force_terms.static"

void stout_smear(stout_control *control, gauge_field_t in)
{
  control->U[0] = in; /* Shallow copy for convenience */
  
  if (!control->calculate_force_terms)
    smear_plain(control, in);
  else
    smear_with_force_terms(control, in);
}
