#include "hex.ih"

void stout_exclude_none(su3_tuple *buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in)
{
  static su3 tmp;

  #define _MULTIPLY_AND_EXPONENTIATE(x, principal) \
  { \
    _real_times_su3(tmp, coeff, (*staples)[x][principal]); \
    _su3_times_su3d(buff_out[x][principal], tmp, buff_in[x][principal]); \
    project_antiherm(&buff_out[x][principal]); \
    exposu3_in_place(&buff_out[x][principal]); \
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
