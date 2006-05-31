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
}

void update_backward_gauge2() {
  int ix=0, kb=0, kb2=0, iz=0, iy=0;

  for(ix = 0; ix < VOLUME/4; ix++) {
    kb2=g_idn[ g_eo2lexic[2*ix] ][0];
    _su3_assign(g_gauge_field_copy[ix][0],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[2*ix]][1];
    _su3_assign(g_gauge_field_copy[ix][1],g_gauge_field[kb][1]);

    kb2=g_idn[g_eo2lexic[2*ix]][2];
    _su3_assign(g_gauge_field_copy[ix][2],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[2*ix]][3];
    _su3_assign(g_gauge_field_copy[ix][3],g_gauge_field[kb][3]);

    kb2=g_idn[ g_eo2lexic[2*ix+1] ][0];
    _su3_assign(g_gauge_field_copy[ix][4],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[2*ix+1]][1];
    _su3_assign(g_gauge_field_copy[ix][5],g_gauge_field[kb][1]);

    kb2=g_idn[g_eo2lexic[2*ix+1]][2];
    _su3_assign(g_gauge_field_copy[ix][6],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[2*ix+1]][3];
    _su3_assign(g_gauge_field_copy[ix][7],g_gauge_field[kb][3]);
  }
  for(ix = 0; ix < VOLUME/4; ix++) {
    iz = VOLUME/4+ix;
    kb = g_eo2lexic[2*ix];
    _su3_assign(g_gauge_field_copy[iz][0],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[iz][1],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[iz][2],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[iz][3],g_gauge_field[kb][3]);

    kb = g_eo2lexic[2*ix+1];
    _su3_assign(g_gauge_field_copy[iz][4],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[iz][5],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[iz][6],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[iz][7],g_gauge_field[kb][3]);
  }
  for(ix = 0; ix < VOLUME/4; ix++) {
    iz = VOLUME/2+ix;
    iy = (VOLUME+RAND)/2+2*ix;
    kb2=g_idn[ g_eo2lexic[iy] ][0];
    _su3_assign(g_gauge_field_copy[iz][0],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[iy]][1];
    _su3_assign(g_gauge_field_copy[iz][1],g_gauge_field[kb][1]);

    kb2=g_idn[g_eo2lexic[iy]][2];
    _su3_assign(g_gauge_field_copy[iz][2],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[iy]][3];
    _su3_assign(g_gauge_field_copy[iz][3],g_gauge_field[kb][3]);

    iy = (VOLUME+RAND)/2+2*ix+1;
    kb2=g_idn[ g_eo2lexic[iy] ][0];
    _su3_assign(g_gauge_field_copy[iz][4],g_gauge_field[kb2][0]);
    kb=g_idn[g_eo2lexic[iy]][1];
    _su3_assign(g_gauge_field_copy[iz][5],g_gauge_field[kb][1]);

    kb2=g_idn[g_eo2lexic[iy]][2];
    _su3_assign(g_gauge_field_copy[iz][6],g_gauge_field[kb2][2]);
    kb=g_idn[g_eo2lexic[iy]][3];
    _su3_assign(g_gauge_field_copy[iz][7],g_gauge_field[kb][3]);
  }

  for(ix = 0; ix < VOLUME/4; ix++) {
    iz = VOLUME/2+VOLUME/4+ix;
    iy = (VOLUME+RAND)/2+2*ix;
    kb = g_eo2lexic[iy];
    _su3_assign(g_gauge_field_copy[iz][0],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[iz][1],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[iz][2],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[iz][3],g_gauge_field[kb][3]);

    kb = g_eo2lexic[iy+1];
    _su3_assign(g_gauge_field_copy[iz][4],g_gauge_field[kb][0]);
    _su3_assign(g_gauge_field_copy[iz][5],g_gauge_field[kb][1]);
    _su3_assign(g_gauge_field_copy[iz][6],g_gauge_field[kb][2]);
    _su3_assign(g_gauge_field_copy[iz][7],g_gauge_field[kb][3]);
  }
}
