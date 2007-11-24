
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
#include "solver/matrix_mult_typedef_nd.h"
#include "sub_low_ev.h"
#include "cg_her.h"

/* P output = solution , Q input = source */
int cg_her_nd(spinor * const P_up,spinor * P_dn, spinor * const Q_up, spinor * const Q_dn, const int max_iter, 
	   double eps_sq, const int rel_prec, const int N, matrix_mult_nd f, 
	   const int subtract_ev, const int modulo){
  double normsp, normsq, pro, err, alpha_cg, beta_cg, squarenorm;
  int iteration;
  double err1, err2;
  
  squarenorm = square_norm(Q_up, N);
  squarenorm+= square_norm(Q_dn, N);
  /*        !!!!   INITIALIZATION    !!!! */
  assign(g_chi_up_spinor_field[DUM_SOLVER], P_up, N);
  assign(g_chi_dn_spinor_field[DUM_SOLVER], P_dn, N);
  
  /*        (r_0,r_0)  =  normsq         */
  normsp =square_norm(P_up, N);
  normsp+=square_norm(P_dn, N);

  if((subtract_ev == 1)) { 
    /* assign_sub_lowest_eigenvalues(g_chi__spinor_field[DUM_SOLVER+5], Q, 10, N); */
  }
  else{
    assign(g_chi_up_spinor_field[DUM_SOLVER+5], Q_up, N);
    assign(g_chi_dn_spinor_field[DUM_SOLVER+5], Q_dn, N);
  }
  
  /* initialize residue r and search vector p */
  if(normsp==0){
    /* if a starting solution vector equal to zero is chosen */
    assign(g_chi_up_spinor_field[DUM_SOLVER+1], g_chi_up_spinor_field[DUM_SOLVER+5], N);
    assign(g_chi_dn_spinor_field[DUM_SOLVER+1], g_chi_dn_spinor_field[DUM_SOLVER+5], N);
    assign(g_chi_up_spinor_field[DUM_SOLVER+2], g_chi_up_spinor_field[DUM_SOLVER+5], N);
    assign(g_chi_dn_spinor_field[DUM_SOLVER+2], g_chi_dn_spinor_field[DUM_SOLVER+5], N);
    normsq =square_norm(Q_up, N);
    normsq+=square_norm(Q_dn, N);
  }
  else{
    /* if a starting solution vector different from zero is chosen */
    f(g_chi_up_spinor_field[DUM_SOLVER+3],g_chi_dn_spinor_field[DUM_SOLVER+3],
      g_chi_up_spinor_field[DUM_SOLVER],g_chi_dn_spinor_field[DUM_SOLVER]);
   
    if((subtract_ev == 1)) {
      /* sub_lowest_eigenvalues(g_chi_spinor_field[DUM_SOLVER+3], g_chi_spinor_field[DUM_SOLVER], 10, N); */
    }
    diff(g_chi_up_spinor_field[DUM_SOLVER+1], g_chi_up_spinor_field[DUM_SOLVER+5], g_chi_up_spinor_field[DUM_SOLVER+3], N);
    diff(g_chi_dn_spinor_field[DUM_SOLVER+1], g_chi_dn_spinor_field[DUM_SOLVER+5], g_chi_dn_spinor_field[DUM_SOLVER+3], N);
    assign(g_chi_up_spinor_field[DUM_SOLVER+2], g_chi_up_spinor_field[DUM_SOLVER+1], N);
    assign(g_chi_dn_spinor_field[DUM_SOLVER+2], g_chi_dn_spinor_field[DUM_SOLVER+1], N);
    normsq =square_norm(g_chi_up_spinor_field[DUM_SOLVER+2], N);
    normsq+=square_norm(g_chi_dn_spinor_field[DUM_SOLVER+2], N);
  }


  /* insert "tmp.c" */



  /* main loop */
  for(iteration=0;iteration<max_iter;iteration++){
    f(g_chi_up_spinor_field[DUM_SOLVER+4],g_chi_dn_spinor_field[DUM_SOLVER+4],
      g_chi_up_spinor_field[DUM_SOLVER+2],g_chi_dn_spinor_field[DUM_SOLVER+2]);

    if((subtract_ev == 1) && (iteration%modulo == 0)) {
      /* sub_lowest_eigenvalues(g_chi_spinor_field[DUM_SOLVER+4], g_chi_spinor_field[DUM_SOLVER+2], 10, N); */
    }
    /* c=scalar_prod(&g_ev[0*VOLUME], g_chi_spinor_field[DUM_SOLVER+4]);
       printf("%e, %e\n",c.re,c.im); */
    pro =scalar_prod_r(g_chi_up_spinor_field[DUM_SOLVER+2], g_chi_up_spinor_field[DUM_SOLVER+4], N);
    pro+=scalar_prod_r(g_chi_dn_spinor_field[DUM_SOLVER+2], g_chi_dn_spinor_field[DUM_SOLVER+4], N);
     
    /*  Compute alpha_cg(i+1)   */
    alpha_cg=normsq/pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(g_chi_up_spinor_field[DUM_SOLVER], g_chi_up_spinor_field[DUM_SOLVER+2],  alpha_cg, N);
    assign_add_mul_r(g_chi_dn_spinor_field[DUM_SOLVER], g_chi_dn_spinor_field[DUM_SOLVER+2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(g_chi_up_spinor_field[DUM_SOLVER+1], g_chi_up_spinor_field[DUM_SOLVER+4], -alpha_cg, N);
    assign_add_mul_r(g_chi_dn_spinor_field[DUM_SOLVER+1], g_chi_dn_spinor_field[DUM_SOLVER+4], -alpha_cg, N);

    /* Check whether the precision is reached ... */
    err1 =square_norm(g_chi_up_spinor_field[DUM_SOLVER+1], N);
    err2 =square_norm(g_chi_dn_spinor_field[DUM_SOLVER+1], N);
    err = err1 + err2;
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("cg_her_nd : i = %d  esqr  %e = %e + %e \n",iteration,err, err1, err2); fflush( stdout);
    }

    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      if((subtract_ev == 1)){
	/* assign_add_invert_subtracted_part(g_chi_spinor_field[DUM_SOLVER], Q, 10, N); */
      } 
      assign(P_up, g_chi_up_spinor_field[DUM_SOLVER], N);
      assign(P_dn, g_chi_dn_spinor_field[DUM_SOLVER], N);
      g_sloppy_precision = 0;
      return(iteration+1);
    }
#ifdef _USE_HALFSPINOR
    if(((err*err <= eps_sq) && (rel_prec == 0)) || ((err*err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      g_sloppy_precision = 1;
      if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
	printf("sloppy precision on\n"); fflush( stdout);
      }
    }
#endif
    /* Compute beta_cg(i+1)
       Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
    beta_cg=err/normsq;
    assign_mul_add_r(g_chi_up_spinor_field[DUM_SOLVER+2], beta_cg, g_chi_up_spinor_field[DUM_SOLVER+1], N);
    assign_mul_add_r(g_chi_dn_spinor_field[DUM_SOLVER+2], beta_cg, g_chi_dn_spinor_field[DUM_SOLVER+1], N);
    normsq=err;
  }
  if((subtract_ev == 1)) { 
    /* assign_add_invert_subtracted_part(g_chi_spinor_field[DUM_SOLVER], Q, 10, N);
       assign_add_invert_subtracted_part(g_chi_spinor_field[DUM_SOLVER], Q, 10, N); */
  }
  assign(P_up, g_chi_up_spinor_field[DUM_SOLVER], N);
  assign(P_dn, g_chi_dn_spinor_field[DUM_SOLVER], N);
  g_sloppy_precision = 0;  
  
  return(-1);
}

static char const rcsid[] = "$Id$";








