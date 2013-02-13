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
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "geometry_eo.h"
#include "measure_rectangles.h"


double measure_rectangles(const su3 ** const gf) {
  static double res;
#ifdef MPI
  double ALIGN mres;
#endif

#ifdef OMP
#pragma omp parallel
  {
  int thread_num = omp_get_thread_num();
#endif

  int i, j, k, mu, nu;
  su3 ALIGN pr1, pr2, tmp; 
  const su3 *v = NULL , *w = NULL;
  double ALIGN ac, ks, kc, tr, ts, tt;

  kc = 0.0;
  ks = 0.0;
#ifdef OMP
#pragma omp for
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
  kc=(kc+ks)/3.0;
#ifdef OMP
  g_omp_acc_re[thread_num] = kc;
#else
  res = kc;
#endif

#ifdef OMP
  } /* OpenMP parallel closing brace */
  
  res = 0.0;
  for(int i = 0; i < omp_num_threads; ++i)
    res += g_omp_acc_re[i];
#else
#endif
#ifdef MPI
  MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  res = mres;
#endif

  return res;
}
