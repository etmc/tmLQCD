/* $Id$ */

/****************************************************
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

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif 


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "solver/solver.h"
#include "mr.h"

int mr(spinor * const P, spinor * const Q,
       const int max_iter, const double eps_sq,
       const int rel_prec, const int N, matrix_mult f){
  int i=0;
  double norm_r,beta;
  complex alpha;

  f(g_spinor_field[DUM_SOLVER+2], P);
  diff(g_spinor_field[DUM_SOLVER], Q, g_spinor_field[DUM_SOLVER+2], N);
  norm_r=square_norm(g_spinor_field[DUM_SOLVER], N);
  if(g_proc_id == g_stdio_proc) {
    printf("MR iteration= %d  |res|^2= %e\n", i, norm_r); 
    fflush( stdout );
  }
  while((norm_r > eps_sq) && (i < max_iter)){
    i++;
    f(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER]);
    alpha=scalar_prod(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER], N);
    beta=square_norm(g_spinor_field[DUM_SOLVER+1], N);
    _mult_real(alpha, alpha, 1./beta);
    assign_add_mul(P, g_spinor_field[DUM_SOLVER], alpha, N);
    if(i%50 == 0){
      f(g_spinor_field[DUM_SOLVER+2], P);
    }
    else{
      assign_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], alpha, N);
    }

    diff(g_spinor_field[DUM_SOLVER], Q, g_spinor_field[DUM_SOLVER+2], N);
    norm_r=square_norm(g_spinor_field[DUM_SOLVER], N);
    if(g_proc_id == g_stdio_proc) {
      printf("MR iteration= %d  |res|^2= %g\n", i, norm_r); 
      fflush(stdout);
    }
 }
  
  if(norm_r > eps_sq){
    return(-1);
  }
  return(i);
}
