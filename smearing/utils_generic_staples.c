#include "utils.ih"

void generic_staples(gauge_field_t buff_out, gauge_field_t buff_in)
{
  static su3 tmp;

#define _ADD_STAPLES_TO_COMPONENT(x, to, via) \
  { \
    _su3_times_su3d(tmp, buff_in.field[g_iup[x][via]][to], buff_in.field[g_iup[x][to]][via]); \
    _su3_times_su3_acc(buff_out.field[x][via], buff_in.field[x][via], tmp); \
    _su3_times_su3(tmp, buff_in.field[g_idn[x][via]][to], buff_in.field[g_iup[g_idn[x][via]][to]][via]); \
    _su3d_times_su3_acc(buff_out.field[x][via], buff_in.field[g_idn[x][via]][via], tmp); \
  }
  
  for (int x = 0; x < VOLUME; ++x)
  {
    _ADD_STAPLES_TO_COMPONENT(x, 0, 1);
    _ADD_STAPLES_TO_COMPONENT(x, 0, 2);
    _ADD_STAPLES_TO_COMPONENT(x, 0, 3);
    
    _ADD_STAPLES_TO_COMPONENT(x, 1, 0);
    _ADD_STAPLES_TO_COMPONENT(x, 1, 2);
    _ADD_STAPLES_TO_COMPONENT(x, 1, 3);

    _ADD_STAPLES_TO_COMPONENT(x, 2, 0);
    _ADD_STAPLES_TO_COMPONENT(x, 2, 1);
    _ADD_STAPLES_TO_COMPONENT(x, 2, 3);

    _ADD_STAPLES_TO_COMPONENT(x, 3, 0);
    _ADD_STAPLES_TO_COMPONENT(x, 3, 1);
    _ADD_STAPLES_TO_COMPONENT(x, 3, 2);
  }
      
#undef _ADD_STAPLES_TO_COMPONENT
}
