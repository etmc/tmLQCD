#include "hyp.ih"

void APE_project_exclude_one(gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in)
{
  double const coeff_principal = 1.0 - coeff;
  double const coeff_staples   = coeff / 4.0;

#define _ADD_AND_REUNITARIZE(x, principal, component) \
  { \
    _real_times_su3_plus_real_times_su3(buff_out.field_array[component / 4].field[x][component % 4], coeff_principal, buff_in.field[x][principal], coeff_staples, staples.field_array[component / 4].field[x][component % 4]) \
    reunitarize(buff_out.field_array[component / 4].field[x] + (component % 4)); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_AND_REUNITARIZE(x, I0_0, I1_0_1);
    _ADD_AND_REUNITARIZE(x, I0_0, I1_0_2);
    _ADD_AND_REUNITARIZE(x, I0_0, I1_0_3);

    _ADD_AND_REUNITARIZE(x, I0_1, I1_1_0);
    _ADD_AND_REUNITARIZE(x, I0_1, I1_1_2);
    _ADD_AND_REUNITARIZE(x, I0_1, I1_1_3);

    _ADD_AND_REUNITARIZE(x, I0_2, I1_2_0);
    _ADD_AND_REUNITARIZE(x, I0_2, I1_2_1);
    _ADD_AND_REUNITARIZE(x, I0_2, I1_2_3);

    _ADD_AND_REUNITARIZE(x, I0_3, I1_3_0);
    _ADD_AND_REUNITARIZE(x, I0_3, I1_3_1);
    _ADD_AND_REUNITARIZE(x, I0_3, I1_3_2);
  }

#undef _ADD_AND_REUNITARIZE
}
