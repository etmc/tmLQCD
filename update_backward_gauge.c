/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "global.h"
#include "su3.h"
#include "update_backward_gauge.h"


#if defined _USE_HALFSPINOR
void update_backward_gauge(su3 ** const gf) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix=0, kb=0, iy=0;

#ifdef OMP
#pragma omp for
#endif 
  for(ix = 0; ix < VOLUME/2; ix++) {
    iy = (VOLUME+RAND)/2+ix;
    kb = g_idn[ g_eo2lexic[iy] ][0];
    _su3_assign(g_gauge_field_copy[0][ix][0], gf[kb][0]);
    kb = g_idn[ g_eo2lexic[iy] ][1];
    _su3_assign(g_gauge_field_copy[0][ix][1], gf[kb][1]);
    kb = g_idn[ g_eo2lexic[iy] ][2];
    _su3_assign(g_gauge_field_copy[0][ix][2], gf[kb][2]);
    kb = g_idn[ g_eo2lexic[iy] ][3];
    _su3_assign(g_gauge_field_copy[0][ix][3], gf[kb][3]);

    kb = g_idn[ g_eo2lexic[ix] ][0];
    _su3_assign(g_gauge_field_copy[1][ix][0], gf[kb][0]);
    kb = g_idn[ g_eo2lexic[ix] ][1];
    _su3_assign(g_gauge_field_copy[1][ix][1], gf[kb][1]);
    kb = g_idn[ g_eo2lexic[ix] ][2];
    _su3_assign(g_gauge_field_copy[1][ix][2], gf[kb][2]);
    kb = g_idn[ g_eo2lexic[ix] ][3];
    _su3_assign(g_gauge_field_copy[1][ix][3], gf[kb][3]);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

  g_update_gauge_copy = 0;
  return;
}

#elif _USE_TSPLITPAR 

void update_backward_gauge(su3 ** const gf) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix=0, kb=0, kb2=0;

#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < VOLUME/2;ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copyt[ix][0],gf[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copyt[ix][1],gf[kb][0]);

    _su3_assign(g_gauge_field_copys[ix][0],gf[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copys[ix][1],gf[kb][1]);

    _su3_assign(g_gauge_field_copys[ix][2],gf[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copys[ix][3],gf[kb][2]);

    _su3_assign(g_gauge_field_copys[ix][4],gf[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copys[ix][5],gf[kb][3]);
  }
#ifdef OMP
#pragma omp for
#endif
  for(ix = (VOLUME+RAND)/2; ix < (VOLUME+RAND)/2+VOLUME/2;ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copyt[ix][0],gf[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copyt[ix][1],gf[kb][0]);

    _su3_assign(g_gauge_field_copys[ix][0],gf[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copys[ix][1],gf[kb][1]);

    _su3_assign(g_gauge_field_copys[ix][2],gf[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copys[ix][3],gf[kb][2]);

    _su3_assign(g_gauge_field_copys[ix][4],gf[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copys[ix][5],gf[kb][3]);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

  g_update_gauge_copy = 0;
  return;
}

#else

void update_backward_gauge(su3 ** const gf) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int ix=0, kb=0, kb2=0;

#ifdef OMP
#pragma omp for
#endif
  for(ix = 0; ix < VOLUME/2; ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copy[ix][0],gf[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copy[ix][1],gf[kb][0]);

    _su3_assign(g_gauge_field_copy[ix][2],gf[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copy[ix][3],gf[kb][1]);

    _su3_assign(g_gauge_field_copy[ix][4],gf[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copy[ix][5],gf[kb][2]);

    _su3_assign(g_gauge_field_copy[ix][6],gf[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copy[ix][7],gf[kb][3]);
  }
#ifdef OMP
#pragma omp for
#endif
  for(ix = (VOLUME+RAND)/2; ix < (VOLUME+RAND)/2+VOLUME/2; ix++) {
    kb2=g_eo2lexic[ix];
    _su3_assign(g_gauge_field_copy[ix][0],gf[kb2][0]);
    kb=g_idn[g_eo2lexic[ix]][0];
    _su3_assign(g_gauge_field_copy[ix][1],gf[kb][0]);

    _su3_assign(g_gauge_field_copy[ix][2],gf[kb2][1]);
    kb=g_idn[g_eo2lexic[ix]][1];
    _su3_assign(g_gauge_field_copy[ix][3],gf[kb][1]);

    _su3_assign(g_gauge_field_copy[ix][4],gf[kb2][2]);
    kb=g_idn[g_eo2lexic[ix]][2];
    _su3_assign(g_gauge_field_copy[ix][5],gf[kb][2]);

    _su3_assign(g_gauge_field_copy[ix][6],gf[kb2][3]);
    kb=g_idn[g_eo2lexic[ix]][3];
    _su3_assign(g_gauge_field_copy[ix][7],gf[kb][3]);
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

  g_update_gauge_copy = 0;
  return;
}

#endif
