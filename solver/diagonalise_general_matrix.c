/* $Id$ */

/******************************************************
 *
 * subroutine to diagonalise a complex n times n
 * matrix. Input is a complex matrix in _C_ like
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
#include "complex.h"
#include "linalg/lapack.h"
#include "diagonalise_general_matrix.h"

void diagonalise_general_matrix(int n, complex * A, int lda, complex * vl,  
				complex * evalues ) {

  complex *vr = NULL, *temp = NULL, *work = NULL, dummy;
  double * rwork = NULL, * scale = NULL, abnrm, * rcone = NULL, * rconv = NULL;
  int lwork, info, i, j, ilo, ihi;
  
  rwork = malloc(2*n*sizeof(double));
  vr = malloc(n*n*sizeof(complex));
/*   temp = malloc(n*n*sizeof(complex)); */
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
  lwork = (int)(dummy.re);
  work = malloc(lwork * sizeof(complex));
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

