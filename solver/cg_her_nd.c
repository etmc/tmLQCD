/***********************************************************************
 *
 * Copyright (C) 2005 Thomas Chiarappa
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
#include "gettime.h"
#include "solver/matrix_mult_typedef_nd.h"
#include "sub_low_ev.h"
#include "solver_field.h"
#include "cg_her.h"

/* P output = solution , Q input = source */
int cg_her_nd(spinor * const P_up,spinor * P_dn, spinor * const Q_up, spinor * const Q_dn, 
	      const int max_iter, double eps_sq, const int rel_prec, 
	      const int N, matrix_mult_nd f) {
  double normsp, normsq, pro, err, alpha_cg, beta_cg, squarenorm;
  int iteration;
  double err1, err2;
  double atime, etime, flops;
  spinor ** up_field = NULL;
  spinor ** dn_field = NULL;  
  const int nr_sf = 5;
  /* do we really need so many fields??? */
  init_solver_field(&up_field, VOLUMEPLUSRAND, nr_sf);
  init_solver_field(&dn_field, VOLUMEPLUSRAND, nr_sf);

  atime = gettime();
  squarenorm = square_norm(Q_up, N, 1);
  squarenorm+= square_norm(Q_dn, N, 1);
  /*        !!!!   INITIALIZATION    !!!! */
  assign(up_field[0], P_up, N);
  assign(dn_field[0], P_dn, N);
  
  /*        (r_0,r_0)  =  normsq         */
  normsp =square_norm(P_up, N, 1);
  normsp+=square_norm(P_dn, N, 1);

/*   assign(up_field[5], Q_up, N); */
/*   assign(dn_field[5], Q_dn, N); */
  
  /* initialize residue r and search vector p */
  if(normsp==0){
    /* if a starting solution vector equal to zero is chosen */
    assign(up_field[1], Q_up, N);
    assign(dn_field[1], Q_dn, N);
    assign(up_field[2], Q_up, N);
    assign(dn_field[2], Q_dn, N);
    normsq =square_norm(Q_up, N, 1);
    normsq+=square_norm(Q_dn, N, 1);
  }
  else {
    /* if a starting solution vector different from zero is chosen */
    f(up_field[3],dn_field[3],
      up_field[0],dn_field[0]);
   
    diff(up_field[1], Q_up, up_field[3], N);
    diff(dn_field[1], Q_dn, dn_field[3], N);
    assign(up_field[2], up_field[1], N);
    assign(dn_field[2], dn_field[1], N);
    normsq =square_norm(up_field[2], N, 1);
    normsq+=square_norm(dn_field[2], N, 1);
  }

  /* main loop */
  for(iteration=0;iteration<max_iter;iteration++){
    f(up_field[4],dn_field[4],
      up_field[2],dn_field[2]);

    pro =scalar_prod_r(up_field[2], up_field[4], N, 1);
    pro+=scalar_prod_r(dn_field[2], dn_field[4], N, 1);
     
    /*  Compute alpha_cg(i+1)   */
    alpha_cg=normsq/pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(up_field[0], up_field[2],  alpha_cg, N);
    assign_add_mul_r(dn_field[0], dn_field[2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(up_field[1], up_field[4], -alpha_cg, N);
    assign_add_mul_r(dn_field[1], dn_field[4], -alpha_cg, N);

    /* Check whether the precision is reached ... */
    err1 =square_norm(up_field[1], N, 1);
    err2 =square_norm(dn_field[1], N, 1);
    err = err1 + err2;
    if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
      printf("cg_her_nd : i = %d  esqr  %e = %e + %e \n",iteration,err, err1, err2); fflush( stdout);
    }

    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      break;
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
    assign_mul_add_r(up_field[2], beta_cg, up_field[1], N);
    assign_mul_add_r(dn_field[2], beta_cg, dn_field[1], N);
    normsq=err;
  }

  assign(P_up, up_field[0], N);
  assign(P_dn, dn_field[0], N);
  g_sloppy_precision = 0;  

  etime = gettime();
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# CG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iteration, eps_sq, etime-atime); 
  }
  
  finalize_solver(up_field, nr_sf);
  finalize_solver(dn_field, nr_sf);
  if(iteration > max_iter) return(-1);
  return(iteration);
}









