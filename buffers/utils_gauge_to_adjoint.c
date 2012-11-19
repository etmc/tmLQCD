#include "utils.ih"

void gauge_to_adjoint(adjoint_field_t *out, gauge_field_t const in)
{
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      _trace_lambda((*out)[x][mu], in[x][mu]);
      (*out)[x][mu].d1 *= -0.5;  (*out)[x][mu].d2 *= -0.5;
      (*out)[x][mu].d3 *= -0.5;  (*out)[x][mu].d4 *= -0.5;
      (*out)[x][mu].d5 *= -0.5;  (*out)[x][mu].d6 *= -0.5;
      (*out)[x][mu].d7 *= -0.5;  (*out)[x][mu].d8 *= -0.5;
    }
}