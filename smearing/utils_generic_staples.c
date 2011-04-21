#include "utils.ih"

void generic_staples(su3 *buff_out, int x, int mu, su3_tuple *buff_in)
{
  static su3 tmp;

#define _ADD_STAPLES_TO_COMPONENT(to, via) \
  { \
    _su3_times_su3d(tmp, buff_in[g_iup[x][via]][to], buff_in[g_iup[x][to]][via]); \
    _su3_times_su3_acc(*buff_out, buff_in[x][via], tmp); \
    _su3_times_su3(tmp, buff_in[g_idn[x][via]][to], buff_in[g_iup[g_idn[x][via]][to]][via]); \
    _su3d_times_su3_acc(*buff_out, buff_in[g_idn[x][via]][via], tmp); \
  }

  _su3_zero(*buff_out);

  switch (mu)
  {
    case 0:
      _ADD_STAPLES_TO_COMPONENT(0, 1);
      _ADD_STAPLES_TO_COMPONENT(0, 2);
      _ADD_STAPLES_TO_COMPONENT(0, 3);
      break;

    case 1:
      _ADD_STAPLES_TO_COMPONENT(1, 0);
      _ADD_STAPLES_TO_COMPONENT(1, 2);
      _ADD_STAPLES_TO_COMPONENT(1, 3);
      break;

    case 2:
      _ADD_STAPLES_TO_COMPONENT(2, 0);
      _ADD_STAPLES_TO_COMPONENT(2, 1);
      _ADD_STAPLES_TO_COMPONENT(2, 3);
      break;

    case 3:
      _ADD_STAPLES_TO_COMPONENT(3, 0);
      _ADD_STAPLES_TO_COMPONENT(3, 1);
      _ADD_STAPLES_TO_COMPONENT(3, 2);
      break;
  }

#undef _ADD_STAPLES_TO_COMPONENT
}
