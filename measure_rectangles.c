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

/*******************************************************************
 *
 * Here the 1x2 rectangles are implemented
 * for renormalization group improved gauge
 * actions like the DBW2 or the Iwasaki
 * gauge action.
 *
 * 1/3 \sum_{\mu\leq\nu;\mu,nu=1}^4 Tr U^{1x2}
 *
 * author: Carsten Urbach
 *         <urbach@physik.fu-berlin.de>
 *
 *******************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "geometry_eo.h"
#include "measure_rectangles.h"


double measure_rectangles(su3 ** const gf) {
/* there is a reduction on ks and kc, so we need to make these OMP shared
 * also, we keep the static keyword because this function relies on the data
 * retention on ga, as can be seen below (it always returns ga) */
  static double ks,kc,ga;

#ifdef OMP
#define static
#endif

#ifdef MPI
  static double gas;
#endif

  ks=0.0; kc=0.0;

#ifdef OMP
#pragma omp parallel
  {
#endif
  int i, j, k, mu, nu;
  static su3 pr1, pr2, tmp; 
  su3 *v = NULL , *w = NULL;
  static double ac, tr, ts, tt;

#ifdef OMP
#undef static
#endif

  if(g_update_rectangle_energy) {
#ifdef OMP
#pragma omp for reduction(+:kc) reduction(+:ks)
#endif
    for (i = 0; i < VOLUME; i++) {
      for (mu = 0; mu < 4; mu++) {
	for (nu = 0; nu < 4; nu++) { 
	  if(nu != mu) {
	    /*
	      ^
	      |
	      ^
	      |
	      ->
	    */
	    j = g_iup[i][mu];
	    k = g_iup[j][nu];
	    v = &gf[i][mu];
	    w = &gf[j][nu];
	    _su3_times_su3(tmp, *v, *w);
	    v = &gf[k][nu];
	    _su3_times_su3(pr1, tmp, *v);
	    /*
	      ->
	      ^
	      |
	      ^
	      |
	    */
	    j = g_iup[i][nu];
	    k = g_iup[j][nu];
	    v = &gf[i][nu];
	    w = &gf[j][nu];
	    _su3_times_su3(tmp, *v, *w);
	    v = &gf[k][mu];
	    _su3_times_su3(pr2, tmp, *v);
	    
	    /* Trace it */
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    /* 	  printf("i mu nu: %d %d %d, ac = %e\n", i, mu, nu, ac); */
	    /* Kahan summation */
	    tr=ac+kc;
	    ts=tr+ks;
	    tt=ts-ks;
	    ks=ts;
	    kc=tr-tt;
	  }
	}
      }
    }

#ifdef OMP
  } /* if g_update_rectangle_energy */
  } /* OpenMP closing brace */

  /* this construct is unfortunately necessary for OpenMP because the closing brace must
   * come before the calculation of ga */
  if(g_update_rectangle_energy) {
#endif

    ga=(kc+ks)/3.0;
  
#ifdef MPI
    MPI_Allreduce(&ga, &gas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ga = gas;
#endif
    g_update_rectangle_energy = 0;
  }
  return ga;
}
