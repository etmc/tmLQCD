#include "hyp.ih"

void hyp_staples_exclude_one(su3_tuple **buff_out, su3_tuple **buff_in)
{
  static su3 tmp;
  for (int idx = 0; idx < 3; ++idx)
    memset(buff_out[idx], 0, sizeof(su3_tuple) * VOLUMEPLUSRAND); /* Brutal but fast zeroing of buffer... */

#define _ADD_STAPLES_TO_COMPONENT(component, to, excl, via, x) \
  { \
    _su3_times_su3d(tmp, buff_in[I2_ ## to ## _ ## excl ## via / 4][g_iup[x][via]][I2_ ## to ## _ ## excl ## via % 4], buff_in[I2_ ## via ## _ ## to ## excl / 4][g_iup[x][to]][I2_ ## via ## _ ## to ## excl % 4]); \
    _su3_times_su3_acc(buff_out[component / 4][x][component % 4], buff_in[I2_ ## via ## _ ## to ## excl / 4][x][I2_ ## via ## _ ## to ## excl % 4], tmp); \
    _su3_times_su3(tmp, buff_in[I2_ ## to ## _ ## excl ## via / 4][g_idn[x][via]][I2_ ## to ## _ ## excl ## via % 4], buff_in[I2_ ## via ## _ ## to ## excl / 4][g_iup[g_idn[x][via]][to]][I2_ ## via ## _ ## to ## excl % 4]); \
    _su3d_times_su3_acc(buff_out[component / 4][x][component % 4], buff_in[I2_ ## via ## _ ## to ## excl / 4][g_idn[x][via]][I2_ ## via ## _ ## to ## excl % 4], tmp); \
  }

  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_STAPLES_TO_COMPONENT(I1_0_1, 0, 1, 2, x);
    _ADD_STAPLES_TO_COMPONENT(I1_0_1, 0, 1, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_0_2, 0, 2, 1, x);
    _ADD_STAPLES_TO_COMPONENT(I1_0_2, 0, 2, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_0_3, 0, 3, 1, x);
    _ADD_STAPLES_TO_COMPONENT(I1_0_3, 0, 3, 2, x);

    _ADD_STAPLES_TO_COMPONENT(I1_1_0, 1, 0, 2, x);
    _ADD_STAPLES_TO_COMPONENT(I1_1_0, 1, 0, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_1_2, 1, 2, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_1_2, 1, 2, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_1_3, 1, 3, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_1_3, 1, 3, 2, x);

    _ADD_STAPLES_TO_COMPONENT(I1_2_0, 2, 0, 1, x);
    _ADD_STAPLES_TO_COMPONENT(I1_2_0, 2, 0, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_2_1, 2, 1, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_2_1, 2, 1, 3, x);
    _ADD_STAPLES_TO_COMPONENT(I1_2_3, 2, 3, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_2_3, 2, 3, 1, x);

    _ADD_STAPLES_TO_COMPONENT(I1_3_0, 3, 0, 1, x);
    _ADD_STAPLES_TO_COMPONENT(I1_3_0, 3, 0, 2, x);
    _ADD_STAPLES_TO_COMPONENT(I1_3_1, 3, 1, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_3_1, 3, 1, 2, x);
    _ADD_STAPLES_TO_COMPONENT(I1_3_2, 3, 2, 0, x);
    _ADD_STAPLES_TO_COMPONENT(I1_3_2, 3, 2, 1, x);
  }

#undef _ADD_STAPLES_TO_COMPONENT
}
