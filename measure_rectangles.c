/* $Id$ */

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "geometry_eo.h"
#include "measure_rectangles.h"


double measure_rectangles() {
  int i, j, k, mu, nu;
  static su3 pr1, pr2, tmp; 
  su3 *v = NULL , *w = NULL;
  static double ga, ac; 
#ifdef MPI
  static double gas;
#endif
  static double ks, kc, tr, ts, tt;
  kc=0.0; ks=0.0;

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
	  v = &g_gauge_field[i][mu];
	  w = &g_gauge_field[j][nu];
	  _su3_times_su3(tmp, *v, *w);
	  v = &g_gauge_field[k][nu];
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
	  v = &g_gauge_field[i][nu];
	  w = &g_gauge_field[j][nu];
	  _su3_times_su3(tmp, *v, *w);
	  v = &g_gauge_field[k][mu];
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
  ga=(kc+ks)/3.0;
#ifdef MPI
  MPI_Allreduce(&ga, &gas, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return gas;
#else
  return ga;
#endif
}
