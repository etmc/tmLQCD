#include "hyp.ih"

void free_hyp_control(hyp_control *control)
{
  /* control->U[0] should always be a shallow copy and not be freed! */
  return_gauge_field(&control->U[1]);
  afree(control->staples[0]);
  afree(control->staples[1]);
  free(control->staples);
  free(control->U);
  free(control);
}
