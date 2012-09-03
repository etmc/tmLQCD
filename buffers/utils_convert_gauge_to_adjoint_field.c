#include "utils.ih"

void convert_gauge_to_adjoint_field(adjoint_field_t dest, gauge_field_t orig)
{
  su3 *op = (su3*)orig;
  su3adj *dp = (su3adj*)dest;
  
  for (int x = 0; x < 4 * VOLUME; ++x)
    _trace_lambda(dp[x], op[x]);
}
