#include "hex.ih"

void stout_exclude_one(gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in)
{
  static su3 tmp;
  
#define _MULTIPLY_AND_EXPONENTIATE(x, principal, component) \
  { \
    _su3_times_su3d(tmp, staples.field_array[component / 4].field[x][component % 4], buff_in.field[x][principal]); \
    project_traceless_antiherm(&tmp); \
    _real_times_su3(buff_out.field_array[component / 4].field[x][component % 4], coeff, tmp); \
    exposu3_in_place(&buff_out.field_array[component / 4].field[x][component % 4]); \
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
