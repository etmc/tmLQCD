#include "hyp.ih"

void free_hyp_control(hyp_control *control)
{
  /* control->U[0] should always be a shallow copy and not be freed! */
  return_gauge_field(&control->U[1]);
  free(control->U);
  free(control);
}
