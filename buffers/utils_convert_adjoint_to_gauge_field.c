#include "utils.ih"

void convert_adjoint_to_gauge_field(gauge_field_t dest, adjoint_field_t orig)
{
  su3adj *op = (su3adj*)orig;
  su3 *dp = (su3*)dest;
  
  for (int x = 0; x < 4 * VOLUME; ++x)
    _make_su3(dp[x], op[x]);
}