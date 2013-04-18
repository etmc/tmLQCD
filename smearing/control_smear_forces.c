#include "control.ih"

void smear_forces(smearing_control_t *control, adjoint_field_t in)
{
  /* General checks can be done here, without the pseudo-dynamic casting code below. */
  if (!control->smearing_performed)
    fatal_error("Need to smear a gauge field before smearing associated forces.", "smear_forces");
  
  switch (control->type)
  {
    case Identity:
      identity_smear_forces((identity_control*)control->type_control, in);
      control->force_result = ((identity_control*)control->type_control)->force_result;
      break;
    case Stout:
      stout_smear_forces((stout_control*)control->type_control, in);
      control->force_result = ((stout_control*)control->type_control)->force_result;
      break;
    case HEX:
      hex_smear_forces((hex_control*)control->type_control, in);
      control->force_result = ((hex_control*)control->type_control)->force_result;
      break;
  }
}
