#include "hyp.ih"

void APE_project_exclude_two(su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in)
{
  double const coeff_principal = 1.0 - coeff;
  double const coeff_staples   = coeff / 2.0;

#define _ADD_AND_REUNITARIZE(x, principal, component) \
  { \
    _real_times_su3_plus_real_times_su3(buff_out[component / 4][x][component % 4], coeff_principal, buff_in[x][principal], coeff_staples, staples[component / 4][x][component % 4]) \
    reunitarize(buff_out[component / 4][x] + (component % 4)); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_AND_REUNITARIZE(x, I0_0, I2_0_12);
    _ADD_AND_REUNITARIZE(x, I0_0, I2_0_23);
    _ADD_AND_REUNITARIZE(x, I0_0, I2_0_13);

    _ADD_AND_REUNITARIZE(x, I0_1, I2_1_02);
    _ADD_AND_REUNITARIZE(x, I0_1, I2_1_03);
    _ADD_AND_REUNITARIZE(x, I0_1, I2_1_23);

    _ADD_AND_REUNITARIZE(x, I0_2, I2_2_01);
    _ADD_AND_REUNITARIZE(x, I0_2, I2_2_03);
    _ADD_AND_REUNITARIZE(x, I0_2, I2_2_13);

    _ADD_AND_REUNITARIZE(x, I0_3, I2_3_01);
    _ADD_AND_REUNITARIZE(x, I0_3, I2_3_02);
    _ADD_AND_REUNITARIZE(x, I0_3, I2_3_12);
  }

#undef _ADD_AND_REUNITARIZE
}
