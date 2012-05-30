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
#include <complex.h>
#include "linalg_eo.h"
#include "linalg/blas.h"
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

void IteratedClassicalGS_bi(_Complex double v[], double *vnrm, int n, int m, _Complex double A[], 
			 _Complex double work1[], int lda) {
  const double alpha = 0.5;

  double vnrm_old;
  int i, isorth = 0;
  int j;
  _Complex double CMONE, CONE;
  char *fupl_n = "N";

  CMONE = -1.;
  CONE = 1.;

  vnrm_old = sqrt(square_norm_bi((bispinor*) v, n*sizeof(_Complex double)/sizeof(bispinor)));

  for(i = 0; !isorth && i < max_cgs_it_bi; i ++) {

    for(j = 0; j < m; j++){
      work1[j] = scalar_prod_bi((bispinor*)(A+j*lda), (bispinor*) v, n*sizeof(_Complex double)/sizeof(bispinor));
    }
    _FT(zgemv)(fupl_n, &n, &m, &CMONE, A, &lda, work1, &ONE, &CONE, v, &ONE, 1);
    (*vnrm) = sqrt(square_norm_bi((bispinor*) v, n*sizeof(_Complex double)/sizeof(bispinor)));

    isorth=((*vnrm) > alpha*vnrm_old);
    vnrm_old = *vnrm;
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

void ModifiedGS_bi(_Complex double v[], int n, int m, _Complex double A[], int lda){

  int i;
  _Complex double s;

  for (i = 0; i < m; i ++) {
    s = -scalar_prod_bi((bispinor*)(A+i*lda), (bispinor*) v, n*sizeof(_Complex double)/sizeof(bispinor));
    _FT(zaxpy)(&n, &s, A+i*lda, &ONE, v, &ONE); 
  }
}

