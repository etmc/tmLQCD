#include "stout.ih"

void free_stout_control(stout_control *control)
{
  return_gauge_field_array(&control->U);
  return_gauge_field_array(&control->Q);
  if (control->calculate_force_terms)
  {
    return_gauge_field_array(&control->f0);
    return_gauge_field_array(&control->f1);
    return_gauge_field_array(&control->f2);
    return_gauge_field_array(&control->u);
    return_gauge_field_array(&control->v);
  }
  free(control);
}
