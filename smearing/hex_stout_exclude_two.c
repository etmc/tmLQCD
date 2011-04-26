#include "hex.ih"

void stout_exclude_two(su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in)
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
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I2_0_12);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I2_0_23);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_0, I2_0_13);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I2_1_02);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I2_1_03);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_1, I2_1_23);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I2_2_01);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I2_2_03);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_2, I2_2_13);

    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I2_3_01);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I2_3_02);
    _MULTIPLY_AND_EXPONENTIATE(x, I0_3, I2_3_12);
  }

#undef _MULTIPLY_AND_EXPONENTIATE
}
