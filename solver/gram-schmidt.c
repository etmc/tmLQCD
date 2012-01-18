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
#include <math.h>
#include <stdio.h>
#include "su3spinor.h"
#include <complex.h>
#include "linalg_eo.h"
#include "linalg/blas.h"
#ifdef CRAY
#include <fortran.h>
#endif
#include "gram-schmidt.h"

const int max_cgs_it=5;
static int ONE = 1;

/*
 *
 *  Iterated Classical Gram-Schmidt Orthogonalization
 *
 *  Orthogonalizes v with respect to A.
 *
 */

void IteratedClassicalGS_old(_Complex double v[], double *vnrm, int n, int m, _Complex double A[], 
			 _Complex double work1[]) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, n2, isorth = 0;
  _Complex double CMONE, CONE, CZERO;
#ifdef CRAY
  char * cupl_c = "C", *cupl_n = "N";
  _fcd fupl_c, fupl_n;
  fupl_c = _cptofcd(cupl_c, strlen(cupl_c));
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
#else
  char * fupl_c = "C", *fupl_n = "N";
#endif

  n2 = 2*n;
  CMONE = -1.;
  CONE = 1.;
  CZERO = 0.;

#ifdef HAVE_LAPACK
  vnrm_old = _FT(dnrm2)(&n2, (double*) v, &ONE);
  for (i = 0; !isorth && i < max_cgs_it; i ++) {
    _FT(zgemv)(fupl_c, &n, &m, &CONE, A, &n, v, &ONE, &CZERO, work1, &ONE, 1);
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &n, work1, &ONE, &CONE, v, &ONE, 1);

    (*vnrm) = _FT(dnrm2)(&n2, (double*) v, &ONE);

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = *vnrm;
  }
#endif
  if (i >= max_cgs_it) {
/*     errorhandler(400,""); */
  }
}


void IteratedClassicalGS(_Complex double v[], double *vnrm, int n, int m, _Complex double A[], 
			 _Complex double work1[], int lda) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, isorth = 0;
  int j;
  _Complex double CMONE, CONE;
#ifdef CRAY
  char *cupl_n = "N";
  _fcd fupl_n;
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
#else
  char *fupl_n = "N";
#endif

  CMONE = -1.;
  CONE = 1.;

  vnrm_old = sqrt(square_norm((spinor*) v, n*sizeof(_Complex double)/sizeof(spinor), 1));

  for(i = 0; !isorth && i < max_cgs_it; i ++) {

    for(j = 0; j < m; j++){
      work1[j] = scalar_prod((spinor*)(A+j*lda), (spinor*) v, n*sizeof(_Complex double)/sizeof(spinor), 1);
    }
#ifdef HAVE_LAPACK
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &lda, work1, &ONE, &CONE, v, &ONE, 1);
#endif
    (*vnrm) = sqrt(square_norm((spinor*) v, n*sizeof(_Complex double)/sizeof(spinor), 1));

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = *vnrm;
  }
  if (i >= max_cgs_it) {
/*     errorhandler(400,""); */
  }
}

#ifdef WITHLAPH

void IteratedClassicalGS_su3vect(_Complex double v[], double *vnrm, int n, int m, _Complex double A[],
				 _Complex double work1[], int lda) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, isorth = 0;
  int j;
  _Complex double CMONE, CONE;

#ifdef CRAY
  char *cupl_n = "N";
  _fcd fupl_n;
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
#else
  char *fupl_n = "N";
#endif

  CMONE = -1.;
  CONE = 1.;

  vnrm_old = sqrt(square_norm_su3vect((su3_vector*) v, n*sizeof(_Complex double)/sizeof(su3_vector),1));

  for(i = 0; !isorth && i < max_cgs_it; i ++) {

    for(j = 0; j < m; j++){
      work1[j] = scalar_prod_su3vect((su3_vector*)(A+j*lda), (su3_vector*) v, n*sizeof(_Complex double)/sizeof(su3_vector),1);
    }
#ifdef HAVE_LAPACK
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &lda, work1, &ONE, &CONE, v, &ONE, 1);
#endif
    (*vnrm) = sqrt(square_norm_su3vect((su3_vector*) v, n*sizeof(_Complex double)/sizeof(su3_vector),1));

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = *vnrm;
  }
  if (i >= max_cgs_it) {
    /*     errorhandler(400,""); */
  }
}

#endif // WITHLAPH

/*
 *  ModifiedGramSchmidt 
 *
 *  Orthogonlaizes v with respect to span{A[:,1:m]}
 */

void ModifiedGS_old(_Complex double v[], int n, int m, _Complex double A[]){

  int i;
  _Complex double s;

  for (i = 0; i < m; i ++) {
    s = scalar_prod((spinor*)(A+i*n), (spinor*) v, n*sizeof(_Complex double)/sizeof(spinor), 1);
    s = -s;
#ifdef HAVE_LAPACK
    _FT(zaxpy)(&n, &s, A+i*n, &ONE, v, &ONE); 
#endif
  }
}


void ModifiedGS(_Complex double v[], int n, int m, _Complex double A[], int lda) {

  int i;
  _Complex double s;

  for (i = 0; i < m; i ++) {
    s = scalar_prod((spinor*)(A+i*lda), (spinor*) v, n*sizeof(_Complex double)/sizeof(spinor), 1);
    s = -s;
#ifdef HAVE_LAPACK
    _FT(zaxpy)(&n, &s, A+i*lda, &ONE, v, &ONE); 
#endif
  }
}

#ifdef WITHLAPH

void ModifiedGS_su3vect(_Complex double v[], int n, int m, _Complex double A[], int lda) {

  int i;
  _Complex double s;

  for (i = 0; i < m; i ++) {
    s = scalar_prod_su3vect((su3_vector*)(A+i*lda), (su3_vector*) v, n*sizeof(_Complex double)/sizeof(su3_vector),1);
    s = -s;
#ifdef HAVE_LAPACK
    _FT(zaxpy)(&n, &s, A+i*lda, &ONE, v, &ONE);
#endif
  }
}

#endif // WITHLAPH
