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
#include "solver/cg_her.h"
#include "sub_low_ev.h"
#include "solver_field.h"
#include "dfl_projector.h"
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
  /* x_0 = P */
  assign(solver_field[0], P, N);
  normsp = square_norm(P, N, 1);

  /* initialize residue r and search vector p */
  /* r0 = b - A x0 */
  f(solver_field[2], solver_field[0]);
  diff(solver_field[1], Q, solver_field[2], N);
  /* z_0 = M^-1 r_0 */
  // here we could have a preconditioner for Q^2
  assign(solver_field[3], solver_field[1], N);
  /* p_0 = z_0 */
  assign(solver_field[2], solver_field[3], N);
  /* (r_0, z_0) */
  pro2 = scalar_prod_r(solver_field[1], solver_field[3], N, 1);  
  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {
    /* w_i = A p_i */
    f(solver_field[4], solver_field[2]);
    /* (p_i, w_i) */
    pro = scalar_prod_r(solver_field[2], solver_field[4], N, 1);
    /*  Compute alpha_cg   */
    alpha_cg = pro2 / pro;
     
    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(solver_field[0], solver_field[2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) A p_i   */
    assign_add_mul_r(solver_field[1], solver_field[4], -alpha_cg, N);

    /* Check whether the precision is reached ... */
    err = square_norm(solver_field[1], N, 1);
    if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
      printf("PCG %d\t%g\n",iteration,err); fflush( stdout);
    }

    if(((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      assign(P, solver_field[0], N);
      g_sloppy_precision = 0;
      finalize_solver(solver_field, nr_sf);
      return(iteration+1);
    }

    /* z_i+1 = M r_(i+1) */
    // here we could have a preconditioner for Q^2
    //mg_Qsq_precon(solver_field[3], solver_field[1]);
    assign(solver_field[3], solver_field[1], N);
    /* Compute beta_cg(i+1) */
    beta_cg = 1. / pro2;
    // here we might use Polak-Ribiere formula instead of the standard one
    // beta = (z_i+1,r_i+1 - r_i) / (z_i,r_i)
    //pro2 = -alpha_cg*scalar_prod_r(solver_field[4], solver_field[3], N, 1);
    // standard choice
    pro2 = scalar_prod_r(solver_field[1], solver_field[3], N, 1);
    beta_cg *= pro2;
    /* p_(i+1) = z_(i+1) + beta_cg p_i */
    assign_mul_add_r(solver_field[2], beta_cg, solver_field[3], N);
  }
  assign(P, solver_field[0], N);
  finalize_solver(solver_field, nr_sf);
  return(-11);
}









