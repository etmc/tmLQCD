
/**************************************************************************
 *
 * $Id$
 *  
 * File: cg_her.c
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
 * Author: Martin Hasenbusch <Martin.Hasenbusch@desy.de> 2001
 * 
 * adapted by Thomas Chiarappa Feb 2003
 * and Carsten Urbach <urbach@ifh.de> (Projecting out the EV)
 *
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "sub_low_ev.h"
#include "cg_her.h"

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif

/* P output = solution , Q input = source */
int cg_her(spinor * const P, spinor * const Q, const int max_iter, 
       double eps_sq, matrix_mult f, 
       const int subtract_ev, const int modulo){
  double normsp, normsq, pro, err, alpha_cg, beta_cg;
  int iteration;
  int N = VOLUME/2;
  
  /*        !!!!   INITIALIZATION    !!!! */
  assign(spinor_field[DUM_SOLVER], P, N);
  /*        (r_0,r_0)  =  normsq         */
  normsp=square_norm(P, N);

  if((subtract_ev == 1)){ 
    assign_sub_lowest_eigenvalues(spinor_field[DUM_SOLVER+5], Q, 10, N); 
  }
  else{
    assign(spinor_field[DUM_SOLVER+5], Q, N);
  } 
  
  /* initialize residue r and search vector p */
  if(normsp==0){
    /* if a starting solution vector equal to zero is chosen */
    assign(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+5], N);
    assign(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+5], N);
    normsq=square_norm(Q, N);
  }
  else{
    /* if a starting solution vector different from zero is chosen */
    f(spinor_field[DUM_SOLVER+3], spinor_field[DUM_SOLVER]);
   
    if((subtract_ev == 1)){
      sub_lowest_eigenvalues(spinor_field[DUM_SOLVER+3], spinor_field[DUM_SOLVER], 10, N);
    }
    diff(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+5], spinor_field[DUM_SOLVER+3], N);
    assign(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+1], N);
    normsq=square_norm(spinor_field[DUM_SOLVER+2], N);
  }
  
  /* main loop */
  for(iteration=0;iteration<max_iter;iteration++){
    f(spinor_field[DUM_SOLVER+4], spinor_field[DUM_SOLVER+2]);

    if((subtract_ev == 1) && (iteration%modulo == 0)){
      sub_lowest_eigenvalues(spinor_field[DUM_SOLVER+4], spinor_field[DUM_SOLVER+2], 10, N);
    }
    /* c=scalar_prod(&g_ev[0*VOLUME], spinor_field[DUM_SOLVER+4]);
       printf("%e, %e\n",c.re,c.im); */
    pro=scalar_prod_r(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+4], N);
     
    /*  Compute alpha_cg(i+1)   */
    alpha_cg=normsq/pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+4], -alpha_cg, N);

    /* Check whether the precision is reached ... */
    err=square_norm(spinor_field[DUM_SOLVER+1], N);
    _SO(if(g_proc_id == g_stdio_proc){printf("%d\t%g\n",iteration,err); fflush( stdout);});

#ifndef _SOLVER_OUTPUT
/*     if ( (iteration%10) == 0 && g_proc_id == g_stdio_proc ) */
/*       printf("%d\t%g\n",iteration,err); */
#endif
    
    if (err <= eps_sq){
      if((subtract_ev == 1)){
	assign_add_invert_subtracted_part(spinor_field[DUM_SOLVER], Q, 10, N); 
      }

      assign(P, spinor_field[DUM_SOLVER], N);
       
      return(iteration+1);
    }
     
    /* Compute beta_cg(i+1)
       Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
    beta_cg=err/normsq;
    assign_mul_add_r(spinor_field[DUM_SOLVER+2], beta_cg, spinor_field[DUM_SOLVER+1], N);
    normsq=err;
  }
  if((subtract_ev == 1)){
    assign_add_invert_subtracted_part(spinor_field[DUM_SOLVER], Q, 10, N); 
  }
  
  assign(P, spinor_field[DUM_SOLVER], N);  

  return(-1);
}

static char const rcsid[] = "$Id$";








