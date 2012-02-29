#include "utils.ih"

void generic_staples(su3 *out, unsigned int x, unsigned int mu, gauge_field_t in)
{
  static su3 tmp; /* NOTE Look for replacement of static variables somehow? */

  for (unsigned int nu = 0; nu < 4; ++nu)
  {
    if (nu == mu)
      continue;
    
    _su3_times_su3d(tmp, in.field[g_iup[x][nu]][mu], in.field[g_iup[x][mu]][nu]);
    _su3_times_su3_acc(*out, in.field[x][nu], tmp);
    _su3_times_su3(tmp, in.field[g_idn[x][nu]][mu], in.field[g_iup[g_idn[x][nu]][mu]][nu]);
    _su3d_times_su3_acc(*out, in.field[g_idn[x][nu]][nu], tmp);
  }
}
