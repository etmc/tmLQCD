void free_smearing_control(smearing_control *control)
{
  switch (control->type)
  {
   case Identity:
      free_identity_control((identity_control*)(control->control));
      break;
   case Stout:
      free_stout_control((stout_control*)(control->control));
      break;
  }
  free(control); 
}

