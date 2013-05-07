#include "utils.ih"

void generic_staples(su3 *out, unsigned int x, unsigned int mu, gauge_field_t in)
{
  su3 ALIGN tmp;

  memset(out, 0, sizeof(su3));

  for (unsigned int nu = 0; nu < 4; ++nu)
  {
    if (nu == mu)
      continue;
    
    /* Prefetch the indices */
    unsigned int const xpm   = g_iup[x][mu];
    unsigned int const xmn   = g_idn[x][nu];
    unsigned int const xpn   = g_iup[x][nu];
    unsigned int const xpmmn = g_idn[xpm][nu];
    
    _su3_times_su3d(tmp, in[xpn][mu], in[xpm][nu]);
    _su3_times_su3_acc(*out, in[x][nu], tmp);
    _su3_times_su3(tmp, in[xmn][mu], in[xpmmn][nu]);
    _su3d_times_su3_acc(*out, in[xmn][nu], tmp);
  }
}
