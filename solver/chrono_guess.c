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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "linalg_eo.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/lu_solve.h"
#include "solver/chrono_guess.h"


/* N is the number of vectors to be stored maximally */
/* _n is the last added vector                       */
/* index_array holds the indices of all the vectors  */
/*  to avoid copying things around                   */
/* V is the volume                                   */
/* trial is the vector to be added                   */
/* v is the array of vectors                         */

void chrono_add_solution(spinor * const trial, spinor ** const v, int index_array[],
			const int N, int * _n, const int V) {

  double norm = 0.;
  int i;

  if(N > 0) {
    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("# CSG: adding vector %d to the list of length %d\n", (*_n)+1, N);
      fflush(stdout);
    }
    if((*_n) < N) {
      index_array[(*_n)] = *_n;
      (*_n)= (*_n)+1;
      /* normalise vector */
      norm = sqrt(square_norm(trial, V, 1));
      mul_r(v[index_array[(*_n)-1]], 1/norm, trial, V);
    }
    else {
      /* Reorder the index_array */
      /* Keep most recent first  */
      for(i = 1; i < N; i++) {
	index_array[i-1] = index_array[i];
      }
      index_array[N-1] = (index_array[N-2]+1)%N;
      /* and normalise */
      norm = sqrt(square_norm(trial, V, 1));
      mul_r(v[index_array[N-1]], 1/norm, trial, V);
    }
  }

  return;
}

/* index_array, _N, _n, V as explained above         */
/* trial is the guess vector to be returned          */
/* phi is the right hand side of A*x = b to be       */
/*   solved                                          */

int chrono_guess(spinor * const trial, spinor * const phi, spinor ** const v, int index_array[], 
		 const int _N, const int _n, const int V, matrix_mult f) {
  int info = 0;
  int i, j, N=_N, n=_n;
  _Complex double s;
  static int init_csg = 0;
  static _Complex double *bn = NULL;
  static _Complex double *G = NULL;
  int max_N = 20;

  if(N > 0) {
    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("# CSG: preparing  trial vector \n");
      fflush(stdout);
    }
    if(init_csg == 0) {
      init_csg = 1;
      bn = (_Complex double*) malloc(max_N*sizeof(_Complex double));
      G = (_Complex double*) malloc(max_N*max_N*sizeof(_Complex double));
    }

    /* Construct an orthogonal basis */
    for(j = n-1; j > n-2; j--) {
      for(i = j-1; i > -1; i--) {
	s = scalar_prod(v[index_array[j]], v[index_array[i]], V, 1);
	assign_diff_mul(v[index_array[i]], v[index_array[j]], s, V);
	if(g_debug_level > 2) {
	  s = scalar_prod(v[index_array[i]], v[index_array[j]], V, 1);
	  if(g_proc_id == 0) {
	    printf("# CSG: <%d,%d> = %e +i %e \n", i, j, creal(s), cimag(s));fflush(stdout);
	  }
	}
      }
    }
    
    /* Generate "interaction matrix" V^\dagger f V */
    /* We assume that f is hermitian               */
    /* Generate also the right hand side           */
    
    for (j = 0; j < n; j++){
      f(trial, v[index_array[j]]);
      
      /* Only the upper triangular part is stored      */
      for(i = 0; i < j+1; i++){
	G[i*N + j] = scalar_prod(v[index_array[i]], trial, V, 1);  
	if(j != i) {
	  (G[j*N + i]) = conj(G[i*N + j]);
	}
	if(g_proc_id == 0 && g_debug_level > 2) {
	  printf("# CSG: G[%d*N + %d]= %e + i %e  \n", i, j, creal(G[i*N + j]), cimag(G[i*N + j]));
	  fflush(stdout);
	}
      }
      /* The right hand side */
      bn[j] = scalar_prod(v[index_array[j]], phi, V, 1);  
    }
    
    /* Solver G y = bn for y and store it in bn */
    LUSolve(n, G, N, bn);

    /* Construct the new guess vector */
    if(info == 0) {
      mul(trial, bn[n-1], v[index_array[n-1]], V); 
      if(g_proc_id == 0 && g_debug_level > 2) {
	printf("# CSG: bn[%d] = %f %f\n", index_array[n-1], creal(bn[index_array[n-1]]), cimag(bn[index_array[n-1]]));
      }
      for(i = n-2; i > -1; i--) {
	assign_add_mul(trial, v[index_array[i]], bn[i], V);
	if(g_proc_id == 0 && g_debug_level > 2) {
	  printf("# CSG: bn[%d] = %f %f\n", index_array[i], creal(bn[index_array[i]]), cimag(bn[index_array[i]]));
	}
      }
    }
    else {
      assign(trial, phi, V);
    }

    if(g_proc_id == 0 && g_debug_level > 1) {
      printf("# CSG: done! n= %d N=%d \n", n, N);fflush(stdout);
    }
  }
  else {
    if(g_proc_id == 0 && g_debug_level > 2) {
      printf("# CSG: using zero trial vector \n");
      fflush(stdout);
    }
    zero_spinor_field(trial, V);
  }

  return(info);
}
