/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
 *               2003 Thomas Chiarappa
 *               2002,2003,2004,2005,2010 Carsten Urbach
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
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/matrix_mult_typedef.h"
#include "sub_low_ev.h"
#include "poly_precon.h"
#include "solver_field.h"
#include "cg_her.h"

int cg_her(spinor * const P, spinor * const Q, const int max_iter, 
           double eps_sq, const int rel_prec, const int N, matrix_mult f) {

  static double normsq,pro,err,alpha_cg,beta_cg,squarenorm;
  int iteration;
  int save_sloppy = g_sloppy_precision;
  double atime, etime, flops;
  spinor ** solver_field = NULL;
  spinor * stmp;
  const int nr_sf = 3;

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  } 
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf); 
  } 
  /* initialize residue r and search vector p */
  atime = gettime();
  squarenorm = square_norm(Q, N, 1);

  f(solver_field[0], P);  

  diff(solver_field[1], Q, solver_field[0], N);
  assign(solver_field[2], solver_field[1], N);
  normsq=square_norm(solver_field[1], N, 1);

  /* main loop */
  for(iteration = 1; iteration <= max_iter; iteration++) {
    f(solver_field[0], solver_field[2]);
    pro = scalar_prod_r(solver_field[2], solver_field[0], N, 1);
    alpha_cg = normsq / pro;
    assign_add_mul_r(P, solver_field[2], alpha_cg, N);

#if (defined SSE2 || defined SSE3)
    assign_mul_add_r(solver_field[0], -alpha_cg, solver_field[1], N);
    err = square_norm(solver_field[0], N, 1);
#else
    err = assign_mul_add_r_and_square(solver_field[0], -alpha_cg, solver_field[1], N, 1);
#endif

    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("CG: iterations: %d res^2 %e\n", iteration, err);
      fflush(stdout);
    }

    if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      break;
    }
#ifdef _USE_HALFSPINOR
    if(((err*err <= eps_sq) && (rel_prec == 0)) || ((err*err <= eps_sq*squarenorm) && (rel_prec == 1))) {
      g_sloppy_precision = 1;
      if(g_debug_level > 2 && g_proc_id == g_stdio_proc && g_sloppy_precision_flag == 1) {
        printf("sloppy precision on\n"); fflush( stdout);
      }
    }
#endif

    beta_cg = err / normsq;
    assign_mul_add_r(solver_field[2], beta_cg, solver_field[0], N);
    stmp = solver_field[0];
    solver_field[0] = solver_field[1];
    solver_field[1] = stmp;
    normsq = err;
  }
  etime = gettime();
  g_sloppy_precision = save_sloppy;
  /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  /* 2*1608.0 because the linalg is over VOLUME/2 */
  flops = (2*(2*1608.0+2*3*4) + 2*3*4 + iteration*(2.*(2*1608.0+2*3*4) + 10*3*4))*N/1.0e6f;
  if(g_debug_level > 0 && g_proc_id == 0 && N != VOLUME) {
    printf("# CG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iteration, eps_sq, etime-atime); 
    printf("# CG: flopcount (for e/o tmWilson only): t/s: %1.4e mflops_local: %.1f mflops: %.1f\n", 
           etime-atime, flops/(etime-atime), g_nproc*flops/(etime-atime));
  }
  finalize_solver(solver_field, nr_sf);
  if(iteration > max_iter) return(-1);
  return(iteration);
}









