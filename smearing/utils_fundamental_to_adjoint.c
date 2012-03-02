#include "utils.ih"

void fundamental_to_adjoint(adjoint_field_t out, gauge_field_t const in)
{
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
      _trace_lambda(out.field[x][mu], in.field[x][mu]);
}