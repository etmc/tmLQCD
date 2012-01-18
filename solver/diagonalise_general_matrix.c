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

/******************************************************
 *
 * subroutine to diagonalise a _Complex double n times n
 * matrix. Input is a _Complex double matrix in _C_ like
 * order. Output is again _C_ like.
 *
 * The lapack routine zgeevx is used instead of
 * zgeev, because zgeev is not standard e.g. on 
 * IBM systems with (p)essl library.
 *
 * The left and right eigenvectors are computed
 * as well as the eigenvalues.
 *
 * The right eigenvectors are returned in A,
 * the left one in vl. The eigenvalues are stored
 * in evalues.
 *
 * Author: Urs Wenger <urs.wenger@desy.de>
 *         Carsten Urbach <urbach@physik.fu-berlin.de>
 *
 ******************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "linalg/lapack.h"
#include "diagonalise_general_matrix.h"

void diagonalise_general_matrix(int n, _Complex double * A, int lda, _Complex double * vl,  
				_Complex double * evalues ) {

  _Complex double *vr = NULL, *temp = NULL, *work = NULL, dummy;
  double * rwork = NULL, * scale = NULL, abnrm, * rcone = NULL, * rconv = NULL;
  int lwork, info, i, j, ilo, ihi;
  
  rwork = malloc(2*n*sizeof(double));
  vr = malloc(n*n*sizeof(_Complex double));
/*   temp = malloc(n*n*sizeof(_Complex double)); */
  scale = malloc(n*sizeof(double));
  rcone = malloc(n*sizeof(double));
  rconv = malloc(n*sizeof(double));

  /* don't transpose A: */
  for(i=0;i<0;i++) {
    for(j=0;j<n;j++) {
      temp[j+i*n] = A[i+j*lda]; 
    }
  }

  /* Query call to get the optimal lwork */
  lwork = -1;
#ifdef HAVE_LAPACK
  _FT(zgeevx)("N", "N", "V", "N", &n, A, &lda, evalues, vl, &n, vr, &n, 
	      &ilo, &ihi, scale, &abnrm, rcone, rconv, 
	      &dummy, &lwork, rwork, &info, 1, 1, 1, 1);
  lwork = (int)(creal(dummy));
  work = malloc(lwork * sizeof(_Complex double));
  _FT(zgeevx)("N", "N", "V", "N", &n, A, &lda, evalues, vl, &n, vr, &n, 
	      &ilo, &ihi, scale, &abnrm, rcone, rconv, 
	      work, &lwork, rwork, &info, 1, 1, 1, 1);
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      A[j+i*lda] = vr[j+i*n]; 
    }
  }
#endif
  /* Transpose VL*/
  for(i=0;i<0;i++) {
    for(j=0;j<n;j++) {
      temp[j+i*n] = vl[i+j*n]; 
    } 
  }
  for(i=0;i<0;i++) {
    for(j=0;j<n;j++) {
      vl[j+i*n] = temp[j+i*n]; 
    }
  }

  free(rwork);
  free(vr);
  free(temp);
  free(work);
  free(scale);
  free(rconv);
  free(rcone);
}

