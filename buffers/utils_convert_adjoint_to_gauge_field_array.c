#include "utils.ih"

void convert_adjoint_to_gauge_field_array(gauge_field_array_t dest, adjoint_field_array_t orig)
{
  for (int ctr = 0; ctr < orig.length; ++ctr)
    convert_adjoint_to_gauge_field(dest.field_array[ctr], orig.field_array[ctr]);
}