#include "hyp.ih"

void hyp_staples_exclude_two(gauge_field_array_t buff_out, gauge_field_t buff_in)
{
  static su3 tmp;
  for (int idx = 0; idx < 3; ++idx)
    memset(buff_out.field_array[idx].field, 0, sizeof(su3_tuple) * VOLUMEPLUSRAND); /* Brutal but fast zeroing of buffer... */

  #define _ADD_STAPELS_TO_COMPONENT(component, to, via, x) \
  { \
    _su3_times_su3d(tmp, buff_in.field[g_iup[x][via]][to], buff_in.field[g_iup[x][to]][via]); \
    _su3_times_su3_acc(buff_out.field_array[component / 4].field[x][component % 4], buff_in.field[x][via], tmp); \
    _su3_times_su3(tmp, buff_in.field[g_idn[x][via]][to], buff_in.field[g_idn[g_iup[x][to]][via]][via]); \
    _su3d_times_su3_acc(buff_out.field_array[component / 4].field[x][component % 4], buff_in.field[g_idn[x][via]][via], tmp); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_STAPELS_TO_COMPONENT(I2_0_12, 0, 3, x);
    _ADD_STAPELS_TO_COMPONENT(I2_0_13, 0, 2, x);
    _ADD_STAPELS_TO_COMPONENT(I2_0_23, 0, 1, x);

    _ADD_STAPELS_TO_COMPONENT(I2_1_02, 1, 3, x);
    _ADD_STAPELS_TO_COMPONENT(I2_1_03, 1, 2, x);
    _ADD_STAPELS_TO_COMPONENT(I2_1_23, 1, 0, x);

    _ADD_STAPELS_TO_COMPONENT(I2_2_01, 2, 3, x);
    _ADD_STAPELS_TO_COMPONENT(I2_2_03, 2, 1, x);
    _ADD_STAPELS_TO_COMPONENT(I2_2_13, 2, 0, x);

    _ADD_STAPELS_TO_COMPONENT(I2_3_01, 3, 2, x);
    _ADD_STAPELS_TO_COMPONENT(I2_3_02, 3, 1, x);
    _ADD_STAPELS_TO_COMPONENT(I2_3_12, 3, 0, x);
  }

#undef _ADD_STAPELS_TO_COMPONENT
}
