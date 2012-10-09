#include "global.h"

#include "update_backward_gauge.h"
#include "dirty_shameful_business.h"

void ohnohack_remap_g_gauge_field(gauge_field_t gf)
{
  fprintf(stderr, "[DEBUG] START -- ohnohack_remap_g_gauge_field\n"); fflush(stderr);
  g_gauge_field[0] = (su3*)gf[0];
  for(unsigned int i = 1; i < VOLUMEPLUSRAND + g_dbw2rand; ++i)
    g_gauge_field[i] = g_gauge_field[i - 1] + 4;
  update_backward_gauge(gf);
  fprintf(stderr, "[DEBUG] FINISH -- ohnohack_remap_g_gauge_field\n"); fflush(stderr);
}
