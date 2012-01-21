/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
 *               2003 Thomas Chiarappa
 *               2002,2003,2004,2005 Carsten Urbach
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
#include "solver/matrix_mult_typedef.h"
#include "solver_field.h"
#include "solver/mixed_cg_her.h"

/* P output = solution , Q input = source */
int mixed_cg_her(spinor * const P, spinor * const Q, const int max_iter, 
		 double eps_sq, const int rel_prec, const int N, matrix_mult f) {

  int i = 0, iter = 0, j = 0;
  double sqnrm = 0., sqnrm2, squarenorm;
  double pro, err, alpha_cg, beta_cg;
  spinor *x, *delta, *y;
  spinor ** solver_field = NULL;
  const int nr_sf = 6;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
  }
  squarenorm = square_norm(Q, N, 1);
  sqnrm = squarenorm;

  delta = solver_field[3];
  x = solver_field[4];
  y = solver_field[5];
  assign(delta, Q, N);
    
  if(squarenorm > 1.e-7) { 
    /* if a starting solution vector different from zero is chosen */
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm = square_norm(delta, N, 1);
    if(((sqnrm <= eps_sq) && (rel_prec == 0)) || ((sqnrm <= eps_sq*squarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      return(0);
    }
  }

  for(i = 0; i < 20; i++) {

    g_sloppy_precision = 1;
    /* main CG loop in lower precision */
    zero_spinor_field(x, N);
    assign(solver_field[1], delta, N);
    assign(solver_field[2], delta, N);
    sqnrm2 = sqnrm;
    for(j = 0; j <= max_iter; j++) {
      f(solver_field[0], solver_field[2]);
      pro = scalar_prod_r(solver_field[2], solver_field[0], N, 1);
      alpha_cg = sqnrm2 / pro;
      assign_add_mul_r(x, solver_field[2], alpha_cg, N);
    
      assign_mul_add_r(solver_field[0], -alpha_cg, solver_field[1], N);
      err = square_norm(solver_field[0], N, 1);

      if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
	printf("inner CG: %d res^2 %g\n", iter+j, err);
	fflush(stdout);
      }
    
      if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
	break;
      }
      beta_cg = err / sqnrm2;
      assign_mul_add_r(solver_field[2], beta_cg, solver_field[0], N);
      assign(solver_field[1], solver_field[0], N);
      sqnrm2 = err;
    }
    /* end main CG loop */
    iter += j;
    g_sloppy_precision = 0;
    add(P, P, x, N);

    f(y, x);
    diff(delta, delta, y, N);
    sqnrm = square_norm(delta, N, 1);
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("mixed CG: true residue %d\t%g\t\n",iter, sqnrm); fflush( stdout);
    }

    if(((sqnrm <= eps_sq) && (rel_prec == 0)) || ((sqnrm <= eps_sq*squarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      return(iter+i);
    }
    iter++;
  }
  finalize_solver(solver_field, nr_sf);
  return(-1);
}









