/***********************************************************************
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
 *
 **************************************************************************/

/* ************************************************************************
 * Conjugate Gradient for su3 vectors
 * Authors: Luigi Scorzato, Marco Cristoforetti
 * based on cg_her.c
 *
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/matrix_mult_typedef.h"
#include "cg_her_su3vect.h"

#ifdef WITHLAPH

int cg_her_su3vect(su3_vector * const P, su3_vector * const Q, const int max_iter, 
		   double eps_sq, const int rel_prec, const int N,const int tslice,  matrix_mult_su3vect f) {

  static double normsq,pro,err,alpha_cg,beta_cg,squarenorm;
  int iteration;
  int save_sloppy = g_sloppy_precision;
  double atime, etime;


  atime = gettime();
  squarenorm = square_norm_su3vect(Q, N, 1);

  f(g_jacobi_field[0],P,tslice);

  diff_su3vect(g_jacobi_field[1], Q, g_jacobi_field[0], N);
  assign_su3vect(g_jacobi_field[2], g_jacobi_field[1], N);
  normsq=square_norm_su3vect(g_jacobi_field[1], N, 1);
  
  /* main loop */
  for(iteration = 1; iteration <= max_iter; iteration++) {
    f(g_jacobi_field[0], g_jacobi_field[2],tslice);
    pro = scalar_prod_r_su3vect(g_jacobi_field[2], g_jacobi_field[0], N, 1);
    alpha_cg = normsq / pro;
    assign_add_mul_r_su3vect(P, g_jacobi_field[2], alpha_cg, N);
    
    assign_mul_add_r_su3vect(g_jacobi_field[0], -alpha_cg, g_jacobi_field[1], N);
    err=square_norm_su3vect(g_jacobi_field[0], N, 1);

    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("CG: iterations: %d res^2 %e\n", iteration, err);
      fflush(stdout);
    }
    
    if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      break;
    }
    beta_cg = err / normsq;
    assign_mul_add_r_su3vect(g_jacobi_field[2], beta_cg, g_jacobi_field[0], N);
    assign_su3vect(g_jacobi_field[1], g_jacobi_field[0], N);
    normsq = err;
  }
  etime = gettime();
  g_sloppy_precision = save_sloppy;
  /* FLOPS= 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  if(g_debug_level > 0  && g_proc_id == 0) {
    printf("CG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iteration, eps_sq, etime-atime); 
  }
  if(iteration > max_iter) return(-1);
  return(iteration);
}

#endif // WITHLAPH
