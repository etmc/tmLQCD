#include "stout_3d.ih"

void free_stout_3d_control(stout_3d_control *control)
{
  if (!control)
    return;

  return_gauge_field(&control->U[1]);
  free(control->U);
  free(control);
  return;
  
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  /* The scratch fields should be initialized and cleared by the using functions, not elsewhere. */
}
