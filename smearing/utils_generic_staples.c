#include "utils.ih"

void generic_staples(su3 *out, unsigned int x, unsigned int mu, gauge_field_t in)
{
  su3 ALIGN tmp;

  memset(out, 0, sizeof(su3));

  for (unsigned int nu = 0; nu < 4; ++nu)
  {
    if (nu == mu)
      continue;
    _su3_times_su3d(tmp, in[g_iup[x][nu]][mu], in[g_iup[x][mu]][nu]);
    _su3_times_su3_acc(*out, in[x][nu], tmp);
    _su3_times_su3(tmp, in[g_idn[x][nu]][mu], in[g_iup[g_idn[x][nu]][mu]][nu]);
    _su3d_times_su3_acc(*out, in[g_idn[x][nu]][nu], tmp);
  }
}
