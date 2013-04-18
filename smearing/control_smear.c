#include "control.ih"

void smear(smearing_control_t *control, gauge_field_t in)
{
  switch (control->type)
  {
    case Identity:
      identity_smear((identity_control*)control->type_control, in);
      control->result = ((identity_control*)control->type_control)->result;
      break;
    case APE:
      ape_smear((ape_control*)control->type_control, in);
      control->result = ((ape_control*)control->type_control)->result;
      break;
    case HYP:
      hyp_smear((hyp_control*)control->type_control, in);
      control->result = ((hyp_control*)control->type_control)->result;
      break;
    case Stout:
      stout_smear((stout_control*)control->type_control, in);
      control->result = ((stout_control*)control->type_control)->result;
      break;
    case HEX:
      hex_smear((hex_control*)control->type_control, in);
      control->result = ((hex_control*)control->type_control)->result;
      break;
  }
  control->smearing_performed = 1;
}
