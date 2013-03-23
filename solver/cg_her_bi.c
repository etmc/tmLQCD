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
 *  
 * File: cg_her_bi.c
 *
 * CG solver for hermitian f only!
 *
 * The externally accessible functions are
 *
 *
 *   int cg_bi(bispinor * const P, bispinor * const Q, double m, const int subtract_ev)
 *     CG solver for bispinor structure
 *
 *
 *
 *   !!!!!  SO FAR NOT IMPLEMENTED FOR EW-SUBTRACTION  !!!!!!
 *
 *
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
#include "solver/matrix_mult_typedef.h"
#include "sub_low_ev.h"
#include "cg_her_bi.h"
#include "solver_field.h"
#include"solver/matrix_mult_typedef_bi.h"


/* P output = solution , Q input = source */
int cg_her_bi(bispinor * const P, bispinor * const Q, const int max_iter, 
       double eps_sq, const int rel_prec, const int N, matrix_mult_bi f) {
  
  double normsp, normsq, pro, err, alpha_cg, beta_cg, squarenorm;
  int iteration;
  bispinor ** bisolver_field = NULL;
  const int nr_sf = 6;
  
  if(N == VOLUME) {
    init_bisolver_field(&bisolver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_bisolver_field(&bisolver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  squarenorm = square_norm((spinor*)Q, 2*N, 1);  
  /*        !!!!   INITIALIZATION    !!!! */
  assign((spinor*)bisolver_field[0], (spinor*)P, 2*N);
  /*        (r_0,r_0)  =  normsq         */
  normsp=square_norm((spinor*)P, 2*N, 1);
  assign((spinor*)bisolver_field[5], (spinor*)Q, 2*N);
  
  /* initialize residue r and search vector p */
  if(normsp == 0) {
    /* if a starting solution vector equal to zero is chosen */
    assign((spinor*)bisolver_field[1], (spinor*)bisolver_field[5], 2*N);
    assign((spinor*)bisolver_field[2], (spinor*)bisolver_field[5], 2*N);
    normsq=square_norm((spinor*)Q, 2*N, 1);
  }
  else {
    /* if a starting solution vector different from zero is chosen */
    f(bisolver_field[3], bisolver_field[0]);
    diff((spinor*)bisolver_field[1], (spinor*)bisolver_field[5], 
	 (spinor*)bisolver_field[3], 2*N);
    assign((spinor*)bisolver_field[2], (spinor*)bisolver_field[1], 2*N);
    normsq=square_norm((spinor*)bisolver_field[2], 2*N, 1);
  }
  
  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    f(bisolver_field[4], bisolver_field[2]);
    pro=scalar_prod_r((spinor*)bisolver_field[2], (spinor*)bisolver_field[4], 2*N, 1);
     
    /*  Compute alpha_cg(i+1)   */
    alpha_cg=normsq/pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r((spinor*)bisolver_field[0], (spinor*)bisolver_field[2],  alpha_cg, 2*N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r((spinor*)bisolver_field[1], (spinor*)bisolver_field[4], -alpha_cg, 2*N);

    /* Check whether the precision is reached ... */
    err=square_norm((spinor*)bisolver_field[1], 2*N, 1);

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("%d\t%g\n",iteration,err); fflush( stdout);
    }
    
    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      assign((spinor*)P, (spinor*)bisolver_field[0], 2*N);
      finalize_bisolver(bisolver_field, nr_sf);
      return(iteration+1);
    }
     
    /* Compute beta_cg(i+1)
       Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
    beta_cg=err/normsq;
    assign_mul_add_r((spinor*)bisolver_field[2], beta_cg, (spinor*)bisolver_field[1], 2*N);
    normsq=err;
  }
  
  assign((spinor*)P, (spinor*)bisolver_field[0], 2*N);  
  finalize_bisolver(bisolver_field, nr_sf);
  return(-1);
}
