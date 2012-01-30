#include "utils.ih"

void convert_gauge_to_adjoint_field_array(adjoint_field_array_t dest, gauge_field_array_t orig)
{
  for (int ctr = 0; ctr < orig.length; ++ctr)
    convert_gauge_to_adjoint_field(dest.field_array[ctr], orig.field_array[ctr]);
}
