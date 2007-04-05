/* $Id$ */

/**************************************************************************
 *
 *
 *  Iterated Classical Gram-Schmidt Orthogonalization for bispinors
 *
 *  Orthogonalizes v with respect to A.
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 *************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <math.h>
#include <stdio.h>
#include "su3spinor.h"
#include "complex.h"
#include "linalg_eo.h"
#include "linalg/blas.h"
#ifdef CRAY
#include <fortran.h>
#endif
#include "gram-schmidt_bi.h"

const int max_cgs_it_bi=5;
static int ONE = 1;

/*
 *
 *  Iterated Classical Gram-Schmidt Orthogonalization
 *
 *  Orthogonalizes v with respect to A.
 *
 */

void IteratedClassicalGS_bi_old(complex v[], double *vnrm, int n, int m, complex A[], 
			 complex work1[]) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, n2, isorth = 0;
  complex CMONE, CONE, CZERO;
#ifdef CRAY
  char * cupl_c = "C", *cupl_n = "N";
  _fcd fupl_c, fupl_n;
  fupl_c = _cptofcd(cupl_c, strlen(cupl_c));
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
#else
  char * fupl_c = "C", *fupl_n = "N";
#endif

  n2 = 2*n;
  CMONE.re = -1.; CMONE.im=0.;
  CONE.re = 1.; CONE.im=0.;
  CZERO.re = 0.; CZERO.im=0.;

  vnrm_old = _FT(dnrm2)(&n2, (double*) v, &ONE);

  for (i = 0; !isorth && i < max_cgs_it_bi; i ++) {
    _FT(zgemv)(fupl_c, &n, &m, &CONE, A, &n, v, &ONE, &CZERO, work1, &ONE, 1);
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &n, work1, &ONE, &CONE, v, &ONE, 1);

    (*vnrm) = _FT(dnrm2)(&n2, (double*) v, &ONE);

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = (*vnrm);
  }
  if (i >= max_cgs_it_bi) {
/*     errorhandler(400,""); */
  }
}


void IteratedClassicalGS_bi(complex v[], double *vnrm, int n, int m, complex A[], 
			 complex work1[], int lda) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, n2, isorth = 0;
  int j;
  complex CMONE, CONE, CZERO;
#ifdef CRAY
  char *cupl_n = "N";
  _fcd fupl_n;
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
#else
  char *fupl_n = "N";
#endif

  n2 = 2*n;
  CMONE.re = -1.; CMONE.im=0.;
  CONE.re = 1.; CONE.im=0.;
  CZERO.re = 0.; CZERO.im=0.;

  vnrm_old = sqrt(square_norm_bi((bispinor*) v, n*sizeof(complex)/sizeof(bispinor)));

  for(i = 0; !isorth && i < max_cgs_it_bi; i ++) {

    for(j = 0; j < m; j++){
      work1[j] = scalar_prod_bi((bispinor*) (A+j*lda), (bispinor*) v, n*sizeof(complex)/sizeof(bispinor));
    }
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &lda, work1, &ONE, &CONE, v, &ONE, 1);
    (*vnrm) = sqrt(square_norm_bi((bispinor*) v, n*sizeof(complex)/sizeof(bispinor)));

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = (*vnrm);
  }
  if (i >= max_cgs_it_bi) {
/*     errorhandler(400,""); */
  }
}


/*
 *  ModifiedGramSchmidt 
 *
 *  Orthogonlaizes v with respect to span{A[:,1:m]}
 */

void ModifiedGS_bi_old(complex v[], int n, int m, complex A[]){

  int i;
  complex s;

  for (i = 0; i < m; i ++) {
    s = scalar_prod_bi((bispinor*) (A+i*n), (bispinor*) v, n*sizeof(complex)/sizeof(bispinor));
    s.re = -s.re; s.im = -s.im;
    _FT(zaxpy)(&n, &s, A+i*n, &ONE, v, &ONE); 
  }
}


void ModifiedGS_bi(complex v[], int n, int m, complex A[], int lda){

  int i;
  complex s;

  for (i = 0; i < m; i ++) {
    s = scalar_prod_bi((bispinor*) (A+i*lda), (bispinor*) v, n*sizeof(complex)/sizeof(bispinor));
    s.re = -s.re; s.im = -s.im;
    _FT(zaxpy)(&n, &s, A+i*lda, &ONE, v, &ONE); 
  }
}

