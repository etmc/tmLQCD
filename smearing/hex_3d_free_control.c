#include "hex_3d.ih"

void free_hex_3d_control(hex_3d_control *control)
{
  /* Result and U[0] are always going to be a shallow copies and should not be cleared! */
  
  return_gauge_field(&control->U[1]);
  free(control->U);
    
  afree(control->V);
    
  free(control);
  return;
}
