/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "global.h"
#include "su3.h"
#include "update_backward_gauge.h"


void update_backward_gauge() {
  int ix=0, kb=0, kb2=0;

#ifdef _NEW_GEOMETRY
  /* set the backward gauge field */
  for(ix = 0; ix < VOLUME;ix++) {
    kb=g_idn[ix][0];
    _su3_assign(g_gauge_field_back[ix][0],g_gauge_field[kb][0]);
    kb=g_idn[ix][1];
    _su3_assign(g_gauge_field_back[ix][1],g_gauge_field[kb][1]);
    kb=g_idn[ix][2];
    _su3_assign(g_gauge_field_back[ix][2],g_gauge_field[kb][2]);
    kb=g_idn[ix][3];
    _su3_assign(g_gauge_field_back[ix][3],g_gauge_field[kb][3]);
  }
#else
  for(ix = 0; ix < VOLUME/2;ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copy[ix][0],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copy[ix][1],g_gauge_field[kb][0]);

    _su3_assign(g_gauge_field_copy[ix][2],g_gauge_field[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copy[ix][3],g_gauge_field[kb][1]);

    _su3_assign(g_gauge_field_copy[ix][4],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copy[ix][5],g_gauge_field[kb][2]);

    _su3_assign(g_gauge_field_copy[ix][6],g_gauge_field[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copy[ix][7],g_gauge_field[kb][3]);
  }
  for(ix = (VOLUME+RAND)/2; ix < (VOLUME+RAND)/2+VOLUME/2;ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copy[ix][0],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copy[ix][1],g_gauge_field[kb][0]);

    _su3_assign(g_gauge_field_copy[ix][2],g_gauge_field[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copy[ix][3],g_gauge_field[kb][1]);

    _su3_assign(g_gauge_field_copy[ix][4],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copy[ix][5],g_gauge_field[kb][2]);

    _su3_assign(g_gauge_field_copy[ix][6],g_gauge_field[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copy[ix][7],g_gauge_field[kb][3]);
  }
#endif
}
