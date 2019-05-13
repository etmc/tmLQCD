/***********************************************************************
 *
 * Copyright (C) 2004 Andrea Shindler
 *               2013 Carsten Urbach
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
 * Author: Andrea Shindler <shindler@ifh.de> Jan 2004
 * 
 * This is a Multi-Shift CG solver
 * it expects that the shifts fulfil
 *
 * shift[0] < shift[1] < shift{2] < ... < shift[no_shifts-1]
 *
 * in modulus. The code will use shift[i]^2, which are all >0
 *
 * parameters:
 * shifts are given to the solver in solver_params->shifts
 * number of shifts is in solver_params->no_shifts
 * the operator to invert in solver_params->M_ndpsi
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/solver.h"
#include "solver_field.h"
#include "cg_mms_tm.h"
#include <io/params.h>

static spinor * ps_qmms;
static spinor ** ps_mms_solver;
static double * sigma;
static double * zitam1, * zita;
static double * alphas, * betas;

extern int index_start;

static void init_mms_tm(const unsigned int nr, const unsigned int N);
static void free_mms_tm();

/* P output = solution , Q input = source */
int cg_mms_tm(spinor ** const P, spinor * const Q,
		 solver_params_t * solver_params) {

  static double normsq, pro, err, squarenorm;
  int iteration, N = solver_params->sdim, no_shifts = solver_params->no_shifts;
  static double gamma, alpham1;
  spinor ** solver_field = NULL;
  double atime, etime;
  const int nr_sf = 3;

  atime = gettime();
  if(solver_params->sdim == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
    init_mms_tm(no_shifts, VOLUMEPLUSRAND);
  } 
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf); 
    init_mms_tm(no_shifts, VOLUMEPLUSRAND/2);
  } 

  zero_spinor_field(P[0], N);
  alphas[0] = 1.0;
  betas[0] = 0.0;
  sigma[0] = solver_params->shifts[0]*solver_params->shifts[0];
  if(g_proc_id == 0 && g_debug_level > 1) printf("# CGMMS: shift %d is %e\n", 0, sigma[0]);

  for(int im = 1; im < no_shifts; im++) {
    sigma[im] = solver_params->shifts[im]*solver_params->shifts[im] - sigma[0];
    if(g_proc_id == 0 && g_debug_level > 1) printf("# CGMMS: shift %d is %e\n", im, sigma[im]);
    // these will be the result spinor fields
    zero_spinor_field(P[im], N);
    // these are intermediate fields
    assign(ps_mms_solver[im-1], Q, N);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }

  /* currently only implemented for P=0 */
  squarenorm = square_norm(Q, N, 1);
  /* if a starting solution vector equal to zero is chosen */
  assign(solver_field[0], Q, N);
  assign(solver_field[1], Q, N);
  normsq = squarenorm;

  /* main loop */
  for(iteration = 0; iteration < solver_params->max_iter; iteration++) {

    /*   Q^2*p and then (p,Q^2*p)  */
    solver_params->M_psi(solver_field[2], solver_field[1]);
    // add the zero's shift
    assign_add_mul_r(solver_field[2], solver_field[1], sigma[0], N);
    pro = scalar_prod_r(solver_field[1], solver_field[2], N, 1);

    /* For the update of the coeff. of the shifted pol. we need alphas[0](i-1) and alpha_cg(i).
       This is the reason why we need this double definition of alpha */
    alpham1 = alphas[0];

    /* Compute alphas[0](i+1) */
    alphas[0] = normsq/pro;
    for(int im = 1; im < no_shifts; im++) {

      /* Now gamma is a temp variable that corresponds to zita(i+1) */ 
      gamma = zita[im]*alpham1/(alphas[0]*betas[0]*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alphas[0]));

      // Now zita(i-1) is put equal to the old zita(i)
      zitam1[im] = zita[im];
      // Now zita(i+1) is updated 
      zita[im] = gamma;
      // Update of alphas(i) = alphas[0](i)*zita(i+1)/zita(i) 
      alphas[im] = alphas[0]*zita[im]/zitam1[im];

      // Compute xs(i+1) = xs(i) + alphas(i)*ps(i) 
      assign_add_mul_r(P[im], ps_mms_solver[im-1], alphas[im], N); 
      // in the CG the corrections are decreasing with the iteration number increasing
      // therefore, we can remove shifts when the norm of the correction vector
      // falls below a threshold
      // this is useful for computing time and needed, because otherwise
      // zita might get smaller than DOUBLE_EPS and, hence, zero
      if(iteration > 0 && (iteration % 10 == 0) && (im == no_shifts-1)) {
	double sn = square_norm(ps_mms_solver[no_shifts-2], N, 1);
        err = alphas[no_shifts-1]*alphas[no_shifts-1]*sn;
        // while because more than one shift could be converged
	while(((err <= solver_params->squared_solver_prec) && (solver_params->rel_prec == 0)) ||
              ((err <= solver_params->squared_solver_prec*squarenorm) && (solver_params->rel_prec > 0))) {
	  no_shifts--;
	  if(g_debug_level > 2 && g_proc_id == 0) {
	    printf("# CGMMS: at iteration %d removed one shift, %d remaining\n", iteration, no_shifts);
      	  }
          if(no_shifts>1) {
            sn = square_norm(ps_mms_solver[no_shifts-2], N, 1);
            err = alphas[no_shifts-1]*alphas[no_shifts-1]*sn;
          } else {
            break;
          }
	}
      }
    }
    
    /*  Compute x_(i+1) = x_i + alphas[0](i+1) p_i    */
    assign_add_mul_r(P[0], solver_field[1],  alphas[0], N);
    /*  Compute r_(i+1) = r_i - alphas[0](i+1) Qp_i   */
    assign_add_mul_r(solver_field[0], solver_field[2], -alphas[0], N);

    /* Check whether the precision eps_sq is reached */

    err = square_norm(solver_field[0], N, 1);

    if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
      printf("# CGMMS iteration: %d residue: %g\n", iteration, err); fflush( stdout );
    }

    if( ((err <= solver_params->squared_solver_prec) && (solver_params->rel_prec == 0)) ||
        ((err <= solver_params->squared_solver_prec*squarenorm) && (solver_params->rel_prec > 0)) ||
        (iteration == solver_params->max_iter -1) ) {
        break;
    }

    /* Compute betas[0](i+1) = (r(i+1),r(i+1))/(r(i),r(i))
       Compute p(i+1) = r(i+1) + beta(i+1)*p(i)  */
    betas[0] = err/normsq;
    assign_mul_add_r(solver_field[1], betas[0], solver_field[0], N);
    normsq = err;

    /* Compute betas(i+1) = betas[0](i+1)*(zita(i+1)*alphas(i))/(zita(i)*alphas[0](i))
       Compute ps(i+1) = zita(i+1)*r(i+1) + betas(i+1)*ps(i)  */
    for(int im = 1; im < no_shifts; im++) {
      betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
      assign_mul_add_mul_r(ps_mms_solver[im-1], solver_field[0], betas[im], zita[im], N);
    }
  }
  etime = gettime();
  g_sloppy_precision = 0;
  if(iteration == solver_params->max_iter -1) iteration = -1;
  else iteration++;
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# CGMMS (%d shifts): iter: %d eps_sq: %1.4e %1.4e t/s\n", solver_params->no_shifts, iteration, solver_params->squared_solver_prec, etime - atime); 
  }

  finalize_solver(solver_field, nr_sf);
  return(iteration);
}


static unsigned int ini_mms = 0;
static unsigned int mms_nr_allocated = 0;

static void init_mms_tm(const unsigned int _nr, const unsigned int N) {
  if(ini_mms == 0 || _nr > mms_nr_allocated) {
    if(mms_nr_allocated != 0) {
      free_mms_tm();
    }

    sigma = (double*)calloc((_nr), sizeof(double));
    zitam1 = (double*)calloc((_nr), sizeof(double));
    zita = (double*)calloc((_nr), sizeof(double));
    alphas = (double*)calloc((_nr), sizeof(double));
    betas = (double*)calloc((_nr), sizeof(double));

    ps_qmms = (spinor*)calloc(N*_nr,sizeof(spinor));
    ps_mms_solver = (spinor**)calloc(_nr,sizeof(spinor*));

    for(int i = 0; i < _nr; i++) {
      ps_mms_solver[i]=(spinor*)(((unsigned long int)(ps_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*N;
    }
    mms_nr_allocated = _nr;
    ini_mms = 1;
  }
}

static void free_mms_tm() {
  free(sigma);
  free(zitam1);
  free(zita);
  free(alphas);
  free(betas);
  free(ps_qmms);
  free(ps_mms_solver);
  mms_nr_allocated = 0;
  ini_mms = 0;
  return;
}
