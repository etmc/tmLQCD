/***********************************************************************
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
 ***********************************************************************/


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
#include "solver_field.h"
#include "pcg_her.h"

/* P output = solution , Q input = source */
int pcg_her(spinor * const P, spinor * const Q, const int max_iter, 
	    double eps_sq, const int rel_prec, const int N, matrix_mult f) {
  double normsp, pro, pro2, err, alpha_cg, beta_cg, squarenorm;
  int iteration;
  spinor ** solver_field = NULL;
  const int nr_sf = 5;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  squarenorm = square_norm(Q, N, 1);
  /*        !!!!   INITIALIZATION    !!!! */
  assign(solver_field[0], P, N);
  /*        (r_0,r_0)  =  normsq         */
  normsp = square_norm(P, N, 1);

  assign(solver_field[3], Q, N);
  /* initialize residue r and search vector p */
  if(normsp==0){
    /* if a starting solution vector equal to zero is chosen */
    /* r0 */
    assign(solver_field[1], solver_field[3], N);
    /* p0 */
  }
  else{
    /* if a starting solution vector different from zero is chosen */
    /* r0 = b - A x0 */
    f(solver_field[2], solver_field[0]);
    diff(solver_field[1], solver_field[3], solver_field[2], N);
  }
  /* z0 = M^-1 r0 */
  invert_eigenvalue_part(solver_field[3], solver_field[1], 10, N);
  /* p0 = z0 */
  assign(solver_field[2], solver_field[3], N);

  /* Is this really real? */
  pro2 = scalar_prod_r(solver_field[1], solver_field[3], N, 1);  
  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    /* A p */
    f(solver_field[4], solver_field[2]);

    pro = scalar_prod_r(solver_field[2], solver_field[4], N, 1);
    /*  Compute alpha_cg(i+1)   */
    alpha_cg=pro2/pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(solver_field[0], solver_field[2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(solver_field[1], solver_field[4], -alpha_cg, N);

    /* Check whether the precision is reached ... */
    err=square_norm(solver_field[1], N, 1);
    if(g_debug_level > 1 && g_proc_id == g_stdio_proc) {
      printf("%d\t%g\n",iteration,err); fflush( stdout);
    }

    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      assign(P, solver_field[0], N);
      g_sloppy_precision = 0;
      finalize_solver(solver_field, nr_sf);
      return(iteration+1);
    }
#ifdef _USE_HALFSPINOR
    if(((err*err <= eps_sq) && (rel_prec == 0)) || ((err*err <= eps_sq*squarenorm) && (rel_prec == 1)) || iteration > 1400) {
      g_sloppy_precision = 1;
      if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
	printf("sloppy precision on\n"); fflush( stdout);
      }
    }
#endif
    /* z_j */
    beta_cg = 1/pro2;
/*     invert_eigenvalue_part(solver_field[3], solver_field[1], 10, N); */
    /* Compute beta_cg(i+1)
       Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
    pro2 = scalar_prod_r(solver_field[1], solver_field[3], N, 1);
    beta_cg *= pro2;
    assign_mul_add_r(solver_field[2], beta_cg, solver_field[3], N);
  }
  assign(P, solver_field[0], N);
  g_sloppy_precision = 0;
/*   return(-1); */
  finalize_solver(solver_field, nr_sf);
  return(1);
}









