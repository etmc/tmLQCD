#include "global.h"

#include "update_backward_gauge.h"
#include "dirty_shameful_business.h"

void ohnohack_remap_g_gauge_field(gauge_field_t gf)
{
  g_gauge_field[0] = (su3*)gf[0];
  for(unsigned int i = 1; i < VOLUMEPLUSRAND + g_dbw2rand; ++i)
    g_gauge_field[i] = g_gauge_field[i - 1] + 4;
#ifdef G_UPDATE_GAUGE_COPY
  update_backward_gauge(gf);
#endif

}

void ohnohack_remap_df0(adjoint_field_t af)
{
  df0[0] = (su3adj*)af[0];
  for(unsigned int i = 1; i < VOLUMEPLUSRAND + g_dbw2rand; ++i)
    df0[i] = df0[i - 1] + 4;
}
