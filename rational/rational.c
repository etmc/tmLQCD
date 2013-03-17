/***********************************************************************
 *
 * Copyright (C) 2013 Carsten Urbach
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
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "zolotarev.h"
#include "rational.h"

// init a rational approximation in the range [eps:1]
// rat->range[0,1] should be the spectral range of the squared operator
// eps is computed to be range[0]/range[1]
// order is the order n of the rational approximation [n,n]
// ca and cb specify the range of monomials to use (0 to order-1)

int init_rational(rational_t * rat, const unsigned int scale) {
  int order = rat->order;
  double * ars = malloc(2*order*sizeof(double));
  double * ar;
  double pmu, pnu;
  double a = rat->range[0], b = rat->range[1];
  double sb = 1.;
  int ca = rat->crange[0], cb = rat->crange[1];

  // sanity check of input parameters
  if(ca > order-1 || cb > order-1 || ca < 0 || cb < 0 || ca > cb || order < 1||
     b < a || b < 0 || a < 0) {
    fprintf(stderr, "parameters to init_rational out of range\n");
    fprintf(stderr, "ca = %d, cb = %d, order = %d, a = %e, b = %e\n", ca, cb, order, a, b);
    return(-1);
  }
  int np = cb - ca + 1;
  if(scale) {
    sb = sqrt(b);
  }
  rat->np = np;
  if(((rat->mu = (double*)malloc(np*sizeof(double))) == NULL)  ||
     ((rat->rmu = (double*)malloc(np*sizeof(double))) == NULL) ||
     ((rat->nu = (double*)malloc(np*sizeof(double))) == NULL)  ||
     ((rat->rnu = (double*)malloc(np*sizeof(double))) == NULL)) {
    fprintf(stderr, "Could not allocate memory for coefficients in init_rational\n");
    return(-2);
  }
  rat->eps = a/b;

  // compute optimal zolotarev approximation
  zolotarev(order, rat->eps, &rat->A, ars, &rat->delta);
  rat->A /= sb;
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("# rational approximation of order %d generated with max deviation delta = %e\n", rat->order, rat->delta);
  }
  // restrict to relevant coefficients [2*ca:2*cb]
  ar = ars + 2*ca;
  // compute mu[] and nu[] = sqrt(ar), mu: r even, nu: r odd
  for (int i = 0; i < np; i++) {
    rat->mu[np-i-1] = sb*sqrt(ar[2*i + 1]);
    rat->nu[np-i-1] = sb*sqrt(ar[2*i]);
  }
  // compute the partial fraction coefficients rmu and rnu
  for (int i = 0; i < np; i++) {  
    pmu = 1.0;
    pnu = 1.0;

    for (int j = 0; j < np; j++) {
      if (j != i) {
	pmu *= ((ar[2*j]-ar[2*i+1]) / (ar[2*j+1]-ar[2*i+1]));
	pnu *= ((rat->mu[j]-rat->nu[i]) / (rat->nu[j]-rat->nu[i]));
      }
    }

    rat->rmu[np-i-1] = sb*sb*(ar[2*i]-ar[2*i+1])*pmu;
    rat->rnu[i] = (rat->mu[i]-rat->nu[i])*pnu;
  }

  free(ars);
  return(0);
}


int free_rational(rational_t * rat) {
  free(rat->mu);
  free(rat->nu);
  free(rat->rmu);
  free(rat->rnu);
  rat->mu = NULL;
  rat->nu = NULL;
  rat->rmu = NULL;
  rat->rnu = NULL;
  return(0);
}
