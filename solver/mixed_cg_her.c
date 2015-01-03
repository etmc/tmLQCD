/***********************************************************************
 * Copyright (C) 2013 Florian Burger
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
 * 
 *  
 * File: mixed_cg_her.c
 *
 * CG solver for hermitian f only!
 *
 * The externally accessible functions are
 *
 *
 *   int cg(spinor * const P, spinor * const Q, double m, const int subtract_ev)
 *     CG solver
 *
 * input:
 *   m: Mass to be use in D_psi
 *   subtrac_ev: if set to 1, the lowest eigenvectors of Q^2 will
 *               be projected out.
 *   Q: source
 * inout:
 *   P: initial guess and result
 * 
 *
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "operator/tm_operators_32.h"
#include "solver/matrix_mult_typedef.h"
#include "read_input.h"

#include "solver_field.h"
#include "solver/mixed_cg_her.h"
#include "gettime.h"




/* P output = solution , Q input = source */
int mixed_cg_her(spinor * const P, spinor * const Q, const int max_iter, 
		 double eps_sq, const int rel_prec, const int N, matrix_mult f, matrix_mult32 f32) {

  int i = 0, iter = 0, j = 0;
  float sqnrm = 0., sqnrm2, squarenorm;
  float pro, err, alpha_cg, beta_cg;
  double sourcesquarenorm, sqnrm_d, squarenorm_d;
  spinor *delta, *y, *xhigh;
  spinor32 *x, *stmp;
  spinor ** solver_field = NULL;
  spinor32 ** solver_field32 = NULL;  
  const int nr_sf = 3;
  const int nr_sf32 = 4;

  int max_inner_it = mixcg_maxinnersolverit;
  int N_outer = max_iter/max_inner_it;
  //to be on the save side we allow at least 10 outer iterations
  if(N_outer < 10) N_outer = 10;
  
  int save_sloppy = g_sloppy_precision_flag;
  double atime, etime, flops;
  
  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);    
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND, nr_sf32);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND/2, nr_sf32);    
  }


  squarenorm_d = square_norm(Q, N, 1);
  //printf("sqarenorm_d = %f\n", squarenorm_d);
  sourcesquarenorm = squarenorm_d;
  sqnrm_d = squarenorm_d;
  squarenorm_d = square_norm(P, N, 1);
 
  //printf("sqarenorm_d = %f\n", squarenorm_d);
 
  delta = solver_field[0];
  y = solver_field[1];
  xhigh = solver_field[2];
  x = solver_field32[3];   
  assign(delta, Q, N);
  
  /*
  if(g_debug_level > 3){
//printf("In mixed_cg solver...\n");
//printf("volume is: %d\n",N); 
   spinor32 * help_low = solver_field32[0];
   spinor * help_high = solver_field[0];
   assign_to_32(help_low, Q, N);
   assign(help_high, Q, N);
   printf("square_norm(Q_high) = %e\n", square_norm(help_high,N,1));
   printf("square_norm(Q_low) = %e\n", square_norm_32(help_low,N,1));  
   f32(solver_field32[1], help_low);
   f(solver_field[1], help_high);
   
   assign_to_64(solver_field[2], solver_field32[1], N);
   diff(solver_field[3], solver_field[1], solver_field[2], N);
   sqnrm = square_norm(solver_field[3], N, 1);
   printf("Operator 32 test: (square_norm) / (spinor component) = %.8e\n", sqnrm/24.0/VOLUME);
   exit(1);
  }
  */
  
  
  if(squarenorm_d > 1.e-7) { 
    /* if a starting solution vector different from zero is chosen */ 
    printf("We have a non-zero starting solution -> using it\n");
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    printf("sqnrm_d = %f\n", sqnrm_d);
    if(((sqnrm_d <= eps_sq) && (rel_prec == 0)) || ((sqnrm_d <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      finalize_solver32(solver_field32, nr_sf32);      
      return(0);
    }
  }
  
  atime = gettime();
  for(i = 0; i < N_outer; i++) {

    /* main CG loop in lower precision */
    zero_spinor_field_32(x, N);
    zero_spinor_field_32(solver_field32[0], N);   
    assign_to_32(solver_field32[1], delta, N);
    assign_to_32(solver_field32[2], delta, N);
    
    sqnrm = (float) sqnrm_d;
    sqnrm2 = sqnrm;
    
    /*inner CG loop */
    for(j = 0; j <= max_inner_it; j++) {
      
      f32(solver_field32[0], solver_field32[2]); 
      pro = scalar_prod_r_32(solver_field32[2], solver_field32[0], N, 1);
      alpha_cg = sqnrm2 / pro;
      
      assign_add_mul_r_32(x, solver_field32[2], alpha_cg, N);
      
      assign_mul_add_r_32(solver_field32[0], -alpha_cg, solver_field32[1], N);      
      
      err = square_norm_32(solver_field32[0], N, 1);

      if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
	printf("inner CG: %d res^2 %g\n", iter+j, err);
	fflush(stdout);
      }
    
      //if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
      if((err <= mixcg_innereps*sqnrm)|| (j==max_inner_it) ||  ((1.3*err <= eps_sq) && (rel_prec == 0)) || ((1.3*err <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
	break;
      }
      beta_cg = err / sqnrm2;
      assign_mul_add_r_32(solver_field32[2], beta_cg, solver_field32[0], N);
      stmp = solver_field32[0];
      solver_field32[0] = solver_field32[1];
      solver_field32[1] = stmp;
      sqnrm2 = err;
      
      
    }
    /* end inner CG loop */
    iter += j;

    /* we want to apply a true double matrix with f(y,P) -> set sloppy off here*/
    g_sloppy_precision_flag = 0;
    
    /* calculate defect in double precision */
    assign_to_64(xhigh, x, N);    
    add(P, P, xhigh, N);
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    if(g_debug_level > 0 && g_proc_id == 0) {
      printf("mixed CG: last inner residue: %g\t\n", err);
      printf("mixed CG: true residue %d %g\t\n",iter, sqnrm_d); fflush(stdout);
    }
    
    /* here we can reset it to its initial value*/
    g_sloppy_precision_flag = save_sloppy;
    
    if(((sqnrm_d <= eps_sq) && (rel_prec == 0)) || ((sqnrm_d <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
      etime = gettime();     
      /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
      /* 2*1608.0 because the linalg is over VOLUME/2 */
      flops = (2*(2*1608.0+2*3*4) + 2*3*4 + iter*(2.*(2*1608.0+2*3*4) + 10*3*4))*N/1.0e6f;
      if(g_debug_level > 0 && g_proc_id == 0 && N != VOLUME) {
         printf("# mixed CG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iter, eps_sq, etime-atime); 
         printf("# mixed CG: flopcount (for e/o tmWilson only): t/s: %1.4e mflops_local: %.1f mflops: %.1f\n", 
         etime-atime, flops/(etime-atime), g_nproc*flops/(etime-atime));
      }
      finalize_solver(solver_field, nr_sf);
      finalize_solver32(solver_field32, nr_sf32); 
      return(iter+i);
    }
    iter++;
  }
  finalize_solver(solver_field, nr_sf);
  finalize_solver32(solver_field32, nr_sf32); 
  return(-1);
}







