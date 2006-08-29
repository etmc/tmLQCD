/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "global.h"
#include "su3.h"
#include "update_backward_gauge.h"


void update_backward_gauge() {
  int ix=0, kb=0, iy=0;

  for(ix = 0; ix < VOLUME/2; ix++) {
    kb =g_idn[ g_eo2lexic[ix] ][0];
    _su3_assign(g_gauge_field_copy[0][0][ix][0],g_gauge_field[kb][0]);

    kb = g_idn[ g_eo2lexic[ix] ][1];
    _su3_assign(g_gauge_field_copy[0][0][ix][1],g_gauge_field[kb][1]);

    kb = g_idn[ g_eo2lexic[ix] ][2];
    _su3_assign(g_gauge_field_copy[0][0][ix][2],g_gauge_field[kb][2]);

    kb = g_idn[ g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copy[0][0][ix][3],g_gauge_field[kb][3]);

  }
  for(ix = 0; ix < VOLUME/2; ix++) {
    kb = g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copy[0][1][ix][0],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[0][1][ix][1],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[0][1][ix][2],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[0][1][ix][3],g_gauge_field[kb][3]);

  }
  for(ix = 0; ix < VOLUME/2; ix++) {
    iy = (VOLUME+RAND)/2+ix;
    kb = g_idn[ g_eo2lexic[iy] ][0];
    _su3_assign(g_gauge_field_copy[1][0][ix][0],g_gauge_field[kb][0]);

    kb = g_idn[ g_eo2lexic[iy] ][1];
    _su3_assign(g_gauge_field_copy[1][0][ix][1],g_gauge_field[kb][1]);

    kb = g_idn[ g_eo2lexic[iy] ][2];
    _su3_assign(g_gauge_field_copy[1][0][ix][2],g_gauge_field[kb][2]);

    kb=g_idn[g_eo2lexic[iy]][3];
    _su3_assign(g_gauge_field_copy[1][0][ix][3],g_gauge_field[kb][3]);

  }

  for(ix = 0; ix < VOLUME/2; ix++) {
    iy = (VOLUME+RAND)/2+ix;
    kb = g_eo2lexic[iy];
    _su3_assign(g_gauge_field_copy[1][1][ix][0],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[1][1][ix][1],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[1][1][ix][2],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[1][1][ix][3],g_gauge_field[kb][3]);

  }
}
