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
    case APE_3D:
      ape_3d_smear((ape_3d_control*)control->type_control, in);
      control->result = ((ape_3d_control*)control->type_control)->result;
      break;
    case HYP:
      hyp_smear((hyp_control*)control->type_control, in);
      control->result = ((hyp_control*)control->type_control)->result;
      break;
//     case HYP_3D:
//       hyp_3d_smear((hyp_3d_control*)control->type_control, in);
//       control->result = ((hyp_3d_control*)control->type_control)->result;
//       break;
    case Stout:
      stout_smear((stout_control*)control->type_control, in);
      control->result = ((stout_control*)control->type_control)->result;
      break;
    case Stout_3D:
      stout_3d_smear((stout_3d_control*)control->type_control, in);
      control->result = ((stout_3d_control*)control->type_control)->result;
      break;
    case HEX:
      hex_smear((hex_control*)control->type_control, in);
      control->result = ((hex_control*)control->type_control)->result;
      break;
    case HEX_3D:
      hex_3d_smear((hex_3d_control*)control->type_control, in);
      control->result = ((hex_3d_control*)control->type_control)->result;
      break;
    default:
      fatal_error("Smearing type not implemented.", "smear");
  }
  control->smearing_performed = 1;
}
