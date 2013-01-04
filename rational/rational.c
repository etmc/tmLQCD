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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "zolotarev.h"
#include "rational.h"

// apply rational approximation as partial fraction
double apply_R(const int order, const double y, const double A, double * rs, double * as) {
  double x = 1.;
  
  for(int i = 0; i < order; i++) {
    x += rs[i]/(y+as[i]);
  }
  return(A*x);
}

// init a rational approximation of order order in range [a:b]
int init_rational(rational_t * rat, const int order, const double a, const double b) {
  double * ars = malloc(2*order*sizeof(double));
  double pmu, pnu, np = order;

  rat->order = order;
  rat->mu = (double*)malloc(order*sizeof(double));
  rat->rmu = (double*)malloc(order*sizeof(double));
  rat->nu = (double*)malloc(order*sizeof(double));
  rat->rnu = (double*)malloc(order*sizeof(double));
  rat->epssq = a*a/b/b;
  rat->range[0] = a;
  rat->range[1] = b;

  zolotarev(order, rat->epssq, &rat->A, ars, &rat->delta);

  // compute mu[] and nu[] = M*sqrt(ar), mu: r even, nu: r odd
  for (int i = 0; i < np; i++) {
    rat->mu[i] = b * sqrt(ars[2*i + 1]);
    rat->nu[i] = b * sqrt(ars[2*i]);
  }
  // compute the partial fraction coefficients rmu and rnu
  for (int i = 0; i < np; i++) {  
    pmu=1.0;
    pnu=1.0;

    for (int j = 0; j < np; j++) {
      if (j!=i) {
	pmu*=((ars[2*j]-ars[2*i+1])/(ars[2*j+1]-ars[2*i+1]));
	pnu*=((rat->mu[j]-rat->nu[i])/(rat->nu[j]-rat->nu[i]));
      }
    }

    rat->rmu[i]=b*b*(ars[2*i]-ars[2*i+1])*pmu;
    rat->rnu[i]=(rat->mu[i]-rat->nu[i])*pnu;
  }

  free(ars);
  return(0);
}

int main() {
  int order = 10;
  double A, eps = 1.e-4, delta;
  double ra = eps, rb = 1.;
  rational_t rat;
  double * ar = malloc(order*sizeof(double));
  
  init_rational(&rat, order, ra, rb);

  for(int i = 0; i < order; i++) {
    ar[i] = (rat.mu[i]/rb)*(rat.mu[i]/rb);
  }
  for(double y = eps*eps; y < 1.; y += eps) {
    double x = apply_R(rat.order, y, rat.A, rat.rmu, ar);
    printf("%e %e %e\n", y, x, 1./sqrt(y));
  }

  return(0);
}
