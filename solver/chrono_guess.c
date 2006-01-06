/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "su3.h"
#include "start.h"
#include "linalg_eo.h"
#include "linalg/fortran.h"
#include "linalg/lapack.h"
#include "linalg/blas.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/chrono_guess.h"

#ifdef HAVE_LAPACK
static int ONE = 1;
static int MONE = -1;
#endif

void chrono_add_solution(spinor * const trial, spinor ** const v, int index_array[],
			const int N, int * _n, const int V) {

  double norm = 0.;
  int i;

  if(N > 0) {
    if((*_n) < N) {
      index_array[(*_n)] = (*_n);
      (*_n)= (*_n)+1;
      norm = sqrt(square_norm(trial, V));
      mul_r(v[index_array[(*_n)-1]], 1/norm, trial, V);
    }
    else {
      /* Reorder the index_array */
      /* Keep most recent first  */
      for(i = 1; i < N; i++) {
	index_array[i-1] = index_array[i];
      }
      index_array[N-1] = (index_array[N-2]+1)%N;
      norm = sqrt(square_norm(trial, V));
      mul_r(v[index_array[N-1]], 1/norm, trial, V);
    }
  }

  return;
}

int chrono_guess(spinor * const trial, spinor * const phi, spinor ** const v, int index_array[], 
		 const int _N, const int _n, const int V, matrix_mult f) {
  int info = 0;
  int i, j, N=_N, n=_n;
  complex s;
  static int init_csg = 0;
  static complex *bn = NULL;
  static complex *work = NULL;
  static complex *G = NULL;
  static int * ipiv = NULL;
  static int lwork = 0;
  int max_N = 0;

  if(N > 0) {
    if(init_csg == 0) {
      init_csg = 1;
      for(i = 0; i < 4; i++) {
	if(g_csg_N[2*i] > max_N) max_N = g_csg_N[2*i];
      }
#ifdef HAVE_LAPACK
      lwork = (_FT(ilaenv)(&ONE, "zhetrf", "U", &max_N, &MONE, &MONE, &MONE, 6, 1)) * N;
#endif
      work = (complex*) malloc(lwork*sizeof(complex));
      bn = (complex*) malloc(max_N*sizeof(complex));
      G = (complex*) malloc(max_N*max_N*sizeof(complex));
      ipiv = (int*) malloc(max_N*sizeof(int));
    }

    /* Construct an orthogonal basis */
    for(j = n-1; j > n-2; j--) {
      for(i = j-1; i > -1; i--) {
	s = scalar_prod(v[index_array[j]], v[index_array[i]], V);
	assign_diff_mul(v[index_array[i]], v[index_array[j]], s, V);
	if(g_debug_level > 2) {
	  s = scalar_prod(v[index_array[i]], v[index_array[j]], V);
	  if(g_proc_id == 0) {
	    printf("CSG: <%d,%d> = %e +i %e \n", i, j, s.re, s.im);fflush(stdout);
	  }
	}
      }
    }
    
    /* Generate "interaction matrix" V^\dagger f V */
    /* We assume that f is hermitian               */
    /* Generate also the right hand side           */
    
    for (j = 0; j < n; j++){
      f(trial, v[index_array[j]]);
      
      /* Only the upper triangel part is stored      */
      for(i = 0; i < j+1; i++){
	G[j*N + i] = scalar_prod(v[index_array[i]], trial, V);  
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("CSG: G[%d*N + %d]= %e + i %e  \n", j, i, G[j*N + i].re, G[j*N + i].im);fflush(stdout);
	}
      } 
      bn[j] = scalar_prod(v[index_array[j]], phi, V);  
    }
    
#ifdef HAVE_LAPACK
    _FT(zhetrf)("U", &n, G, &N, ipiv, work, &lwork, &info, 1);
#endif
    if(info != 0) {
      printf("Error in zhetrf info = %d\n", info);
    }
    else {
#ifdef HAVE_LAPACK
      _FT(zhetrs)("U", &n, &ONE, G, &N, ipiv, bn, &N, &info, 1);
#endif
      if(info != 0) {
	printf("Error in zhetrs info = %d\n", info);
      }
    }
    
    if(info == 0) {
      mul(trial, bn[n-1], v[index_array[n-1]], V); 
      if(g_proc_id == 0 && g_debug_level > 1) {
	printf("CSG: bn[%d] = %f %f\n", index_array[n-1], bn[index_array[n-1]].re, bn[index_array[n-1]].im);
      }
      for(i = n-2; i > -1; i--) {
	assign_add_mul(trial, v[index_array[i]], bn[i], V);
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("CSG: bn[%d] = %f %f\n", index_array[i], bn[index_array[i]].re, bn[index_array[i]].im);
	}
      }
    }
    else {
      assign(trial, phi, V);
    }

    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("CSG: done! n= %d N=%d \n", n, N);fflush(stdout);
    }
  }
  else {

#if !defined HAVE_LAPACK
    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("CSG: No lapack available -> No CSG \n");fflush(stdout);
    }
#endif
    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("Using zero trial vector \n");
      fflush(stdout);
    }
    zero_spinor_field(trial, V);
  }

  return(info);
}
