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
  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);    
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND, nr_sf32);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND/2, nr_sf32);    
  }

  squarenorm_d = square_norm(Q, N, 1);
  sourcesquarenorm = squarenorm_d;
  sqnrm_d = squarenorm_d;
  
  //maybe that is wrong
  squarenorm_d = square_norm(P, N, 1);
  
  delta = solver_field[0];
  y = solver_field[1];
  xhigh = solver_field[2];
  x = solver_field32[3];  
  assign(delta, Q, N);
  
  /* here comes a small test*/
//   spinor32 * help_low = solver_field32[0];
//   spinor * help_high = solver_field[0];
//   assign_to_32(help_low, Q, N);
//   assign(help_high, Q, N);
//   printf("square_norm(Q_high) = %e\n", square_norm(help_high,N,1));
//   printf("square_norm(Q_low) = %e\n", square_norm_32(help_low,N,1));  
//   f32(solver_field32[1], help_low);
//   f(solver_field[1], help_high);
//   
//   assign_to_64(solver_field[2], solver_field32[1], N);
//   diff(solver_field[3], solver_field[1], solver_field[2], N);
//   sqnrm = square_norm(solver_field[3], N, 1);
//   printf("Operator 32 test: (square_norm) / (spinor component) = %.8e\n", sqnrm/24.0/VOLUME);
//   exit(1);
  /* end of small test*/
  
  
  /* small performance test */
  
//   int Nhit = 40;
//   double t1,t2,dt,sdt;
//   double antioptaway=0.0;
//   t1 = gettime();
//   antioptaway=0.0;
//   for (i=0;i<Nhit;i++) {
//      f(y, delta);
//      f(delta,y);
//   }
//      antioptaway+=creal(delta[0].s0.c0);  
//   t2 = gettime();
//   dt = t2-t1;
//   sdt=2*1.0e6f*dt/((double)(Nhit*(N)));
//   printf("antioptaway = %e\n", antioptaway);
//   printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt), (int)sizeof(spinor)/3);
//   assign(delta, Q, N);
//   
//   
//   assign_to_32(x, delta, N);
//   float antioptaway_f=0.0;
//   t1 = gettime();
//   antioptaway_f=0.0;
//   for (i=0;i<Nhit;i++) {
//      f32(solver_field32[0], x);
//      f32(x,solver_field32[0]);
// 
//   }
//   antioptaway_f+=creal(x[0].s0.c0);
//   t2 = gettime();
//   dt = t2-t1;
//   sdt=2*1.0e6f*dt/((double)(Nhit*(N)));
//   printf("antioptaway = %e\n", antioptaway_f);
//   printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt), (int)sizeof(spinor32)/3);
//   
  /* end of small performance test */
  
    
  if(squarenorm_d > 1.e-7) { 
    /* if a starting solution vector different from zero is chosen */
    printf("We have a non-zero starting solution -> using it\n");
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    if(((sqnrm <= eps_sq) && (rel_prec == 0)) || ((sqnrm <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      finalize_solver32(solver_field32, nr_sf32);      
      return(0);
    }
  }

  
  for(i = 0; i < 20; i++) {

    //g_sloppy_precision = 1;
    //g_sloppy_precision_flag = 1;
    /* main CG loop in lower precision */
    zero_spinor_field_32(x, N);
    zero_spinor_field_32(solver_field32[0], N);   
    assign_to_32(solver_field32[1], delta, N);
    assign_to_32(solver_field32[2], delta, N);
    
    sqnrm = (float) sqnrm_d;
    sqnrm2 = sqnrm;
    for(j = 0; j <= max_iter; j++) {
      f32(solver_field32[0], solver_field32[2]); 
      
      pro = scalar_prod_r_32(solver_field32[2], solver_field32[0], N, 1);
      alpha_cg = sqnrm2 / pro;
      
      assign_add_mul_r_32(x, solver_field32[2], alpha_cg, N);
      
      assign_mul_add_r_32(solver_field32[0], -alpha_cg, solver_field32[1], N);      
      
      err = square_norm_32(solver_field32[0], N, 1);

      if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
	printf("inner CG: %d res^2 %g\n", iter+j, err);
	fflush(stdout);
      }
    
      //if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
      if((err <= 6.0e-6*sqnrm)) {
	break;
      }
      beta_cg = err / sqnrm2;
      assign_mul_add_r_32(solver_field32[2], beta_cg, solver_field32[0], N);
      stmp = solver_field32[0];
      solver_field32[0] = solver_field32[1];
      solver_field32[1] = stmp;
      sqnrm2 = err;
      
    }
    /* end main CG loop */
    iter += j;
    //g_sloppy_precision = 0;
    //g_sloppy_precision_flag = 0;
    assign_to_64(xhigh, x, N);    
    
//     f(y, xhigh);
//     diff(P, delta, y, N);
//     sqnrm = square_norm(P, N, 1);
//     printf("mixed CG: true residue %d\t%g\t\n",iter, sqnrm);   
//     exit(1);
    
    
    add(P, P, xhigh, N);
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("mixed CG: true residue %d\t%g\t\n",iter, sqnrm_d); fflush(stdout);
    }

    if(((sqnrm_d <= eps_sq) && (rel_prec == 0)) || ((sqnrm_d <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
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







