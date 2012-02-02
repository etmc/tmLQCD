#include "utils.ih"

void cayley_hamilton_exponent(gauge_field_t U, gauge_field_t Q)
{
  complex_field_array_t f0 = get_complex_field_array(4);
  complex_field_array_t f1 = get_complex_field_array(4);
  complex_field_array_t f2 = get_complex_field_array(4);
  real_field_array_t u     = get_real_field_array(4);
  real_field_array_t v     = get_real_field_array(4);
  
  cayley_hamilton_exponent_with_force_terms(U, Q, u, v, f0, f1, f2);
}
