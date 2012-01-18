/***********************************************************************
 *
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
 *
 * Minimal residual solver
 * int mr(spinor * const P, spinor * const Q,
 *	const int max_iter, const double eps_sq,
 *	matrix_mult f){ *
 *
 * returns the number of iterations needed to reach
 * the desired precision. return -1 if the maximal
 * number of iterations was reached.
 *
 * Inout:                                                                      
 *  spinor * P       : guess for the solving spinor                                             
 * Input:                                                                      
 *  spinor * Q       : source spinor
 *  int max_iter     : maximal number of iterations                                 
 *  double eps_sqr   : stopping criterium                                                     
 *  matrix_mult f    : pointer to a function containing 
 *                     the matrix mult for type 
 *                     matrix_mult see 
 *                     matrix_mult_typedef.h
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "start.h"
#include "su3.h"
#include "linalg_eo.h"
#include "solver/solver.h"
#include "solver_field.h"
#include "mr.h"

int mr(spinor * const P, spinor * const Q,
       const int max_iter, const double eps_sq,
       const int rel_prec, const int N, const int parallel, 
       matrix_mult f){
  int i=0;
  double norm_r,beta;
  _Complex double alpha;
  spinor * r;
  spinor ** solver_field = NULL;
  const int nr_sf = 3;
  
  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  r = solver_field[0];
  
  zero_spinor_field(P, N);
  f(solver_field[2], P);
  diff(r, Q, solver_field[2], N);
  norm_r=square_norm(solver_field[0], N, parallel);
  if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
    printf("MR iteration number: %d, |res|^2 = %e\n", i, norm_r); 
    fflush( stdout );
  }
  while((norm_r > eps_sq) && (i < max_iter)){
    i++;
    f(solver_field[1], r);
    alpha=scalar_prod(solver_field[1], r, N, parallel);
    beta=square_norm(solver_field[1], N, parallel);
    alpha /= beta;
    assign_add_mul(P, r, alpha, N);
    if(i%50 == 0){
      f(solver_field[2], P);
    }
    else{
      assign_add_mul(solver_field[2], solver_field[1], alpha, N);
    }

    diff(r, Q, solver_field[2], N);
    norm_r=square_norm(solver_field[0], N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("# MR iteration= %d  |res|^2= %g\n", i, norm_r); 
      fflush(stdout);
    }
  }
  finalize_solver(solver_field, nr_sf);
  if(norm_r > eps_sq){
    return(-1);
  }
  return(i);
}


int mrblk(spinor * const P, spinor * const Q,
	  const int max_iter, const double eps_sq,
	  const int rel_prec, const int N, 
	  matrix_mult_blk f, const int blk) {
  static int mr_init=0;
  int i = 0;
  double norm_r,beta;
  _Complex double alpha;
  spinor * r;
  const int parallel = 0;
  spinor * s[3];
  static spinor *s_=NULL;
  static int N_;

  if(mr_init == 0 || N != N_) {
    if(N!= N_ && mr_init != 0) {
      free(s_);
    }
    N_ = N;
    s_ = calloc(3*(N+1)+1, sizeof(spinor));
    mr_init = 1;
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  s[0] = (spinor *)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE); 
#else
  s[0] = s_;
#endif
  s[1] = s[0] + N + 1;
  s[2] = s[1] + N + 1;

  r = s[0];
  norm_r = square_norm(Q, N, parallel);
  
  zero_spinor_field(P, N);
  f(s[2], P, blk);
  diff(r, Q, s[2], N);
  norm_r = square_norm(r, N, parallel);
  if(g_proc_id == g_stdio_proc && g_debug_level > 2 && blk == 0) {
    printf("MRblk iteration= %d  |res|^2= %e\n", i, norm_r);
    fflush( stdout );
  }
  
  while((norm_r > eps_sq) && (i < max_iter)){
    i++;
    f(s[1], r, blk);
    alpha = scalar_prod(s[1], r, N, parallel);
    beta = square_norm(s[1], N, parallel);
    alpha /= beta;
    assign_add_mul(P, r, alpha, N);
    if(i%50 == 0) {
      f(s[2], P,blk);
    }
    else{
      assign_add_mul(s[2], s[1], alpha, N);
    }
    
    diff(r, Q, s[2], N);
    norm_r = square_norm(r, N, parallel);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2 && blk == 0) {
      printf("MRblk iteration= %d  |res|^2= %g\n", i, norm_r);
      fflush(stdout);
    }
  }
  /* free(s_); */
  if(norm_r > eps_sq){
    return(-1);
  }
  return(i);
}
