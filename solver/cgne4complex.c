/***********************************************************************
 * Copyright (C) 2013 Carsten Urbach
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
# include<tmlqcd_config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"global.h"
#include"su3.h"
#include"linalg_eo.h"
#include"gettime.h"
#include"gcr4complex.h"
#include"cgne4complex.h"

int cgne4complex(_Complex double * const P, _Complex double * const Q, 
		 const int max_iter, const double eps_sq, const int rel_prec,
		 const int N, const int lda, c_matrix_mult f) {
  
  double normsq, pro, err, alpha_cg, beta_cg, squarenorm;
  _Complex double *w_f[3], * _w_f, *stmp;
  double atime, etime;
  int iter;
  
  _w_f = (_Complex double *)malloc(3*lda*sizeof(_Complex double));
  w_f[0] = _w_f; w_f[1] = _w_f+lda; w_f[2] = _w_f+2*lda;
  
    /* initialize residue r and search vector p */
  atime = gettime();
  squarenorm = lsquare_norm(Q, N, 1);

  f(w_f[0], P);  

  ldiff(w_f[1], Q, w_f[0], N);
  memcpy(w_f[2], w_f[1], N*sizeof(_Complex double));
  normsq=lsquare_norm(w_f[1], N, 1);

  /* main loop */
  for(iter = 1; iter <= max_iter; iter++) {
    f(w_f[0], w_f[2]);
    pro = lscalar_prod_r(w_f[2], w_f[0], N, 1);
    alpha_cg = normsq / pro;
    lassign_add_mul_r(P, w_f[2], alpha_cg, N);

    lassign_mul_add_r(w_f[0], -alpha_cg, w_f[1], N);
    err = lsquare_norm(w_f[0], N, 1);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("lCG: iterations: %d res^2 %e\n", iter, err);
      fflush(stdout);
    }

    if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      break;
    }

    beta_cg = err / normsq;
    lassign_mul_add_r(w_f[2], beta_cg, w_f[0], N);
    stmp = w_f[0];
    w_f[0] = w_f[1];
    w_f[1] = stmp;
    normsq = err;
  }
  etime = gettime();
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# lCG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iter, eps_sq, etime-atime); 
  }
  free(_w_f);
  if(iter > max_iter) return(-1);
  return(iter);

}
