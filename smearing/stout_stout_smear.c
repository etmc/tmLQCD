#include "stout.ih"

#include "stout_smear_plain.static"
#include "stout_smear_with_force_terms.static"

void stout_smear(struct stout_control *control, gauge_field_t in)
{
  if (!control->calculate_force_terms)
    stout_smear_plain(control, in);
  else
    stout_smear_with_force_terms(control, in);
}
