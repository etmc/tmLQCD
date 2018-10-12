#include "hex.ih"
#include "global.h"

void stout_exclude_one(su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in)
{
  static su3 tmp;

#define _MULTIPLY_AND_EXPONENTIATE(x, principal, component) \
  { \
    _su3_times_su3d(tmp, staples[component / 4][x][component % 4], buff_in[x][principal]); \
    project_antiherm(&tmp); \
    _real_times_su3(buff_out[component / 4][x][component % 4], coeff, tmp); \
    exposu3_in_place(&buff_out[component / 4][x][component % 4]); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I1_0_1);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I1_0_2);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I1_0_3);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I1_1_0);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I1_1_2);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I1_1_3);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I1_2_0);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I1_2_1);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I1_2_3);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I1_3_0);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I1_3_1);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I1_3_2);
  }

#undef _MULTIPLY_AND_EXPONENTIATE
}
