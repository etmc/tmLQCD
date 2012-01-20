/***********************************************************************
 *
 * Copyright (C) 2007 Carsten Urbach
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
 *
 *
 * Hybrid-Monte-Carlo for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

const double epss=1.e-7;
const int MT = 10;
#define MR 8

inline dmax(double x, double y) {
  if(x > y) return(x);
  return(y);
}

int laguer(double complex a[], const int m, double complex *x, int *its, const int maxit) {
  int iter, i, j;
  double abx, abp, abm, err;
  double complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
  static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
  for (iter = 1; iter <= maxit; iter++) { 
    *its = iter;
    b = a[m];
    err = cabs(b);
    d = 0.;
    f = 0.;
    abx = cabs(*x);
    for (j = m-1; j >= 0; j--) {
      f = (*x) * f + d;
      d = (*x) * d + b;
      b = (*x) * b + a[j];
      err = cabs(b) + abx * err;
    }
    err *= epss;
    if (cabs(b) <= err) return(0);
    g = d / b;
    g2 = g * g;
    h = g2 - 2. * f / b;
    sq = csqrt((double)(m-1) * ((double)(m)*h - g2));
    gp = g + sq;
    gm = g - sq;
    abp = cabs(gp);
    abm = cabs(gm);
    if (abp < abm) gp = gm;
    dx=((dmax(abp,abm) > 0. ?
	 ((double complex)(m))/gp :
	 (1. + abx)*(cos((double)iter) + _Complex_I*sin((double)iter))));
    x1 = (*x) - dx;
    
    if (creal(*x) == creal(x1) && cimag(*x) == cimag(x1)) {
      return(0);
    }
    if (iter % MT) {
      *x=x1;
    }
    else {
      *x = (*x) - frac[iter/MT]*dx;
    }
  }
  fprintf(stderr, "Too many iterations in laguer\n");
  return(-1);
}

int zroots(double complex a[], const int m, double complex roots[], const int polish) {
  int i, j, jj, its, k;
  double complex x, b, c, ad[1000];
  for(j = 0; j < m+1; j++) {
    ad[j] = a[j];
  }
  for(j = m; j > 0; j--) {
    x = 0.;
    if((k = laguer(ad, j, &x, &its, 800)) != 0) {
      fprintf(stderr, "something wront!\n");
    }
    if(abs(cimag(x)) <= 2.*epss*abs(creal(x))) x = creal(x);
    roots[j-1] = x;
    b = ad[j];
    for(jj = j-1; jj > -1; jj--) {
      c = ad[jj];
      ad[jj] = b;
      c = x*b + c;
    }
  }
  if(polish) {
    for(j = 1; j < m+1; j++) {
      if((k = laguer(a, m, &roots[j-1], &its, 800)) != 0) {
	fprintf(stderr, "something wront!\n");
      }
    }
  }
  for(j = 2; j < m+1; j++) {
    x = roots[j-1];
    for(i = j-1; i > 0; i--) {
      if(creal(roots[i-1]) <= creal(x)) break;
      roots[i] = roots[i-1];
    }
  }
  return(0);
}

int main() {
  int i;
  double complex a[5];
  double complex roots[5];

  a[0] = 1.;
  a[1] = 1.;
  a[2] = 1.;
  a[3] = 1.;
  a[4] = 1.;
  zroots(a, 2, roots, 1);
  for(i = 0; i < 2; i++) {
    printf("%f %f %f %f\n", creal(roots[i]), cimag(roots[i]), 
	   creal(a[0]+a[1]*roots[i]+a[2]*roots[i]*roots[i]),
	   cimag(a[0]+a[1]*roots[i]+a[2]*roots[i]*roots[i]));
  }

  return(0);
}
