/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "linalg/fortran.h"
#include "linalg/lapack.h"
#include "linalg/blas.h"
#include "solver/matrix_mult_typedef.h"
/* #include "tm_operators.h" */
#include "solver/chrono_guess.h"

static int ONE = 1;
static int MONE = -1;

int chrono_guess(spinor * const trial, spinor * const phi, spinor ** const v, int index_array[], 
		 const int _N, int * _n, const int V, matrix_mult f) {
  int i, j, n, info, N=_N;
  complex s;
  static int init_csg = 0;
  static complex *bn = NULL;
  static complex *work = NULL;
  static complex *G = NULL;
  static int * ipiv = NULL;
  static int lwork = 0;
  
  if(init_csg == 0) {
    init_csg = 1;
    lwork = (_FT(ilaenv)(&ONE, "zhetrf", "U", &N, &MONE, &MONE, &MONE, 6, 1)) * N;
    work = (complex*) malloc(lwork*sizeof(complex));
    bn = (complex*) malloc(N*sizeof(complex));
    G = (complex*) malloc(N*N*sizeof(complex));
    ipiv = (int*) malloc(N*sizeof(int));
  }


  if((*_n) < N) {
    index_array[(*_n)] = (*_n);
    (*_n)= (*_n)+1;
  }
  else {
    /* Reorder the index_array */
    /* Keep most recent first  */
    for(i = 1; i < N; i++) {
      index_array[i-1] = index_array[i];
    }
  }

  n = (*_n);
  printf("CSG: n= %d N=%d id=%d\n", n, N, g_proc_id);fflush(stdout);

  /* Construct an orthogonal basis */
  for(i = n-2; i > -1; i--) {
    s = scalar_prod(v[index_array[i]], trial, V);
    s.re = -s.re; s.im = -s.im;
    _FT(zaxpy)(&n, &s, (complex*) v[index_array[i]], &ONE, (complex*) trial, &ONE); 
  }
  assign(v[index_array[n-1]], trial, V);

  printf("CSG: orthogonalized n= %d N=%d id=%d\n", n, N, g_proc_id);fflush(stdout);
  
  /* Generate "interaction matrix" V^\dagger f V */
  /* We assume that f is hermitian               */
  /* Generate also the right hand side           */

  for (j = 0; j < n; j++){
    printf("CSG: loop G_mn j= %d index_array[%d] = %d id=%d \n", j, j, index_array[j], g_proc_id);fflush(stdout);
    f(trial, v[index_array[j]]);
    printf("CSG: after f loop G_mn j= %d id = %d\n", j, g_proc_id);fflush(stdout);

    /*     idummy = j+1; */
    /* Only the upper triangel part is stored      */
    for(i = 0; i < j+1; i++){
      if(g_proc_id == 0) {
	printf("CSG: before G loop G_mn j= %d \n", j);fflush(stdout);
      } 
      G[j*N + i] = scalar_prod(v[index_array[i]], trial, V); 
      if(g_proc_id == 0) {
	printf("CSG: loop G_mn j= %d i= %d\n", j, i);fflush(stdout);
      }
    } 
    /*     _FT(zgemv)(fupl_c, &n, &idummy, &CONE, V, &n, temp1, &ONE, */
    /* 	       &CZERO, M+i*jmax, &ONE, 1); */
    bn[j] = scalar_prod(v[index_array[j]], phi, V);
  }

  if(g_proc_id == 0) {
    printf("CSG: computed G_mn n= %d N=%d \n", n, N);fflush(stdout);
  }

  _FT(zhetrf)("U", &n, G, &N, ipiv, work, &lwork, &info, 1);
  if(info != 0) {
    printf("Error in zhetrf info = %d\n", info);
  }
  else {
    _FT(zhetrs)("U", &n, &ONE, G, &N, ipiv, bn, &N, &info, 1);
    if(info != 0) {
      printf("Error in zhetrs info = %d\n", info);
    }
  }
  
  if(g_proc_id == 0) {
    printf("CSG: solved n= %d N=%d \n", n, N);fflush(stdout);
  }

  if(info == 0) {
    mul(trial, bn[index_array[n-1]], v[index_array[n-1]], V);
    for(i = n-2; i > -1; i--) {
      assign_add_mul(trial, v[index_array[i]], bn[index_array[i]], V);
    }
  }
  else {
    assign(trial, phi, V);
  }

  if(g_proc_id == 0) {
    printf("CSG: done! n= %d N=%d \n", n, N);fflush(stdout);
  }
  
  return(0);
}
