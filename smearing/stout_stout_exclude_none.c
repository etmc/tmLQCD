#include "stout.ih"

void stout_exclude_none(gauge_field_t buff_out, double const coeff, gauge_field_t staples, gauge_field_t buff_in)
{
  static su3 tmp;

#define _MULTIPLY_AND_EXPONENTIATE(x, principal) \
  { \
    _su3_times_su3d(tmp, staples.field[x][principal], buff_in.field[x][principal]); \
    project_antiherm(&tmp); \
    _real_times_su3(buff_out.field[x][principal], coeff, tmp); \
    exposu3_in_place(&buff_out.field[x][principal]); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_1);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_2);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_3);
  }

#undef _MULTIPLY_AND_EXPONENTIATE
}
