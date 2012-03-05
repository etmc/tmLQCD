#include "utils.ih"

void adjoint_to_gauge(gauge_field_t out, adjoint_field_t const in)
{
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
      _make_su3(out.field[x][mu], in.field[x][mu]);
}