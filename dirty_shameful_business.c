#include "dirty_shameful_business.h"

void ohnohack_remap_g_gauge_field(su3_tuple *gf)
{
  g_gauge_field[0] = (su3*)gf[0];
  for(i = 1; i < V; i++)
    g_gauge_field[i] = g_gauge_field[i - 1] + 4;
  update_backward_gauge(gf);
}
