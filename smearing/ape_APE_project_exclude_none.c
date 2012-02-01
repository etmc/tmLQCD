#include "ape.ih"

void APE_project_exclude_none(gauge_field_t buff_out, double const coeff, gauge_field_t staples, gauge_field_t buff_in)
{
  double const coeff_principal = 1.0 - coeff;
  double const coeff_staples   = coeff / 6.0;

#define _ADD_AND_REUNITARIZE(x, component) \
  { \
    _real_times_su3_plus_real_times_su3(buff_out.field[x][component], coeff_principal, buff_in.field[x][component], coeff_staples, staples.field[x][component]) \
    reunitarize(buff_out.field[x] + component); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_AND_REUNITARIZE(x, I0_0);
    _ADD_AND_REUNITARIZE(x, I0_1);
    _ADD_AND_REUNITARIZE(x, I0_2);
    _ADD_AND_REUNITARIZE(x, I0_3);
  }

#undef _ADD_AND_REUNITARIZE
}
