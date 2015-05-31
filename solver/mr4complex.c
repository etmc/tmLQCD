/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2010 claude Tadonki
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
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"gcr4complex.h"
#include"mr4complex.h"

int mr4complex(_Complex double * const P, _Complex double * const Q,
	       const int max_iter, const double eps_sq,
	       const int rel_prec, const int N, 
	       const int parallel, const int lda,
	       c_matrix_mult f){
  int i=0;
  double norm_r, norm_sq, beta;
  _Complex double alpha;
  _Complex double *w_f[3], * _w_f, *stmp;
  _Complex double *r;
  
  _w_f = (_Complex double *)malloc(3*lda*sizeof(_Complex double));
  w_f[0] = _w_f; w_f[1] = _w_f+lda; w_f[2] = _w_f+2*lda;
  
  r = w_f[0];
  
  for(int j = 0; j < N; j++) {
    P[j] = 0.;
    w_f[2][j] = 0.;
    r[j] = Q[j];
  }
  //f(w_f[2], P);
  //ldiff(r, Q, w_f[2], N);
  norm_sq = lsquare_norm(Q, N, parallel);
  norm_r = lsquare_norm(r, N, parallel);
  if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
    printf("lMR iteration number: %d, |res|^2 = %e of %e %d\n", i, norm_r, eps_sq, rel_prec); 
    fflush( stdout );
  }
  while(((norm_r > eps_sq && !rel_prec) || ((norm_r > eps_sq*norm_sq && rel_prec))) && (i < max_iter)) {
    i++;
    f(w_f[1], r);
    alpha=lscalar_prod(w_f[1], r, N, parallel);
    beta=lsquare_norm(w_f[1], N, parallel);
    alpha /= beta;
    lassign_add_mul(P, r, alpha, N);
    if(i%50 == 0){
      f(w_f[2], P);
    }
    else{
      lassign_add_mul(w_f[2], w_f[1], alpha, N);
    }

    ldiff(r, Q, w_f[2], N);
    norm_r=lsquare_norm(w_f[0], N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("# lMR iteration= %d  |res|^2= %g\n", i, norm_r); 
      fflush(stdout);
    }
  }
  if(norm_r > eps_sq){
    return(-1);
  }
  return(i);
}
