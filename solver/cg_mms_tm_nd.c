/***********************************************************************
 *
 *
 * Copyright (C) 2004 Andrea Shindler
 *               2009 Carsten Urbach
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
#include "solver/solver.h"
#include "solver_field.h"
#include "cg_mms_tm_nd.h"

static spinor * ps_qmms;
static spinor ** ps_mms_solver;
static double * sigma;
static double * zitam1, * zita;
static double * alphas, * betas;

extern int index_start;

static void init_mms_tm_nd(const unsigned int nr, const unsigned int N);
static void free_mms_tm_nd();

/* P output = solution , Q input = source */
int cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_pm_t * solver_pm) {

  static double normsq, pro, err, squarenorm;
  int iteration, N = solver_pm->N;
  static double gamma, alpham1;
  spinor ** solver_field = NULL;
  const int nr_sf = 4;

  if(solver_pm->N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, 2*nr_sf);
  } 
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, 2*nr_sf); 
  } 

  // don't need boundaries, because we never apply f to them
  // so N is enough
  init_mms_tm_nd(solver_pm->no_shifts, solver_pm->N);

  /* currently only implemented for P=0 */
  for(int im = 0; im < solver_pm->no_shifts; im++) {
    sigma[im] = solver_pm->shifts[im]*solver_pm->shifts[im];
    // these will be the result spinor fields
    zero_spinor_field(Pup[im], N);
    zero_spinor_field(Pdn[im], N);
    // these are intermediate fields
    assign(ps_mms_solver[2*im], Qup, N);
    assign(ps_mms_solver[2*im+1], Qdn, N);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }

  squarenorm = square_norm(Qup, N, 1) + square_norm(Qdn, N, 1);
  /* if a starting solution vector equal to zero is chosen */
  assign(solver_field[0], Qup, N);
  assign(solver_field[1], Qdn, N);
  assign(solver_field[2], Qup, N);
  assign(solver_field[3], Qdn, N);
  normsq = squarenorm;

  /* main loop */
  for(iteration = 0; iteration < solver_pm->max_iter; iteration++) {

    /*   Q^2*p and then (p,Q^2*p)  */
    solver_pm->g(solver_field[6], solver_field[7], solver_field[2], solver_field[3]);
    // add the zero's shift
    assign_add_mul_r(solver_field[6], solver_field[2], sigma[0], N);
    assign_add_mul_r(solver_field[7], solver_field[3], sigma[0], N);
    pro = scalar_prod_r(solver_field[2], solver_field[6], N, 1);
    pro += scalar_prod_r(solver_field[3], solver_field[7], N, 1);

    /* For the update of the coeff. of the shifted pol. we need alphas[0](i-1) and alpha_cg(i).
       This is the reason why we need this double definition of alpha */
    alpham1 = alphas[0];

    /* Compute alphas[0](i+1) */
    alphas[0] = normsq/pro;
    for(int im = 1; im < solver_pm->no_shifts; im++) {

      /* Now gamma is a temp variable that corresponds to zita(i+1) */ 
      gamma = zita[im]*alpham1/(alphas[0]*betas[0]*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alphas[0]));

      /* Now zita(i-1) is put equal to the old zita(i) */
      zitam1[im] = zita[im];
      /* Now zita(i+1) is updated */
      zita[im] = gamma;
      /* Update of alphas(i) = alphas[0](i)*zita(i+1)/zita(i) */ 
      alphas[im] = alphas[0]*zita[im]/zitam1[im];
      /* Compute xs(i+1) = xs(i) + alphas(i)*ps(i) */
      assign_add_mul_r(Pup[im], ps_mms_solver[2*im], alphas[im], N); 
      assign_add_mul_r(Pdn[im], ps_mms_solver[2*im+1], alphas[im], N); 
    }

    /*  Compute x_(i+1) = x_i + alphas[0](i+1) p_i    */
    assign_add_mul_r(Pup[0], solver_field[2],  alphas[0], N);
    assign_add_mul_r(Pdn[0], solver_field[3],  alphas[0], N);
    /*  Compute r_(i+1) = r_i - alphas[0](i+1) Qp_i   */
    assign_add_mul_r(solver_field[0], solver_field[6], -alphas[0], N);
    assign_add_mul_r(solver_field[1], solver_field[7], -alphas[0], N);

    /* Check whether the precision eps_sq is reached */

    err = square_norm(solver_field[0], N, 1) + square_norm(solver_field[1], N, 1);

    if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
      printf("# CGMMS iteration: %d residue: %g\n", iteration, err); fflush( stdout );
    }

    if( ((err <= solver_pm->eps_sq) && (solver_pm->rel_prec == 0)) ||
      ((err <= solver_pm->eps_sq*squarenorm) && (solver_pm->rel_prec == 1)) ) {

      if(g_debug_level > 2) {
	solver_pm->g(solver_field[2], solver_field[3], Pup[0], Pdn[0]);
	diff(solver_field[4], solver_field[2], Qup, N);
	diff(solver_field[5], solver_field[3], Qdn, N);
	err = square_norm(solver_field[4], N, 1) + square_norm(solver_field[5], N, 1);
	if(g_proc_id == g_stdio_proc) {
	  printf("# CG MMS true residue at final iteration (%d) was %g.\n", iteration, err); 
	  fflush( stdout);
	}
      }
      g_sloppy_precision = 0;

      /* save all the results of (Q^dagger Q)^(-1) \gamma_5 \phi */
      /* here ... */
      /* when im == -1 save the base mass*/
      for(int im = 1; im < solver_pm->no_shifts; im++) {
	betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
	assign_mul_add_mul_r(ps_mms_solver[2*im], solver_field[0], betas[im], zita[im], N);
	assign_mul_add_mul_r(ps_mms_solver[2*im+1], solver_field[1], betas[im], zita[im], N);
      }
      finalize_solver(solver_field, nr_sf);
      return(iteration+1);
    }

    /* Compute betas[0](i+1) = (r(i+1),r(i+1))/(r(i),r(i))
       Compute p(i+1) = r(i+1) + beta(i+1)*p(i)  */
    betas[0] = err/normsq;
    assign_mul_add_r(solver_field[2], betas[0], solver_field[0], N);
    assign_mul_add_r(solver_field[3], betas[0], solver_field[1], N);
    normsq = err;

    /* Compute betas(i+1) = betas[0](i)*(zita(i+1)*alphas(i))/(zita(i)*alphas[0](i))
       Compute ps(i+1) = zita(i+1)*r(i+1) + betas(i+1)*ps(i)  */
    for(int im = 1; im < solver_pm->no_shifts; im++) {
      betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
      assign_mul_add_mul_r(ps_mms_solver[2*im], solver_field[0], betas[im], zita[im], N);
      assign_mul_add_mul_r(ps_mms_solver[2*im+1], solver_field[1], betas[im], zita[im], N);
    }
  }
  g_sloppy_precision = 0;
  finalize_solver(solver_field, nr_sf);
  return(-1);
}


static unsigned int ini_mms_nd = 0;
static unsigned int nr_nd = 0;

static void init_mms_tm_nd(const unsigned int _nr, const unsigned int N) {
  if(ini_mms_nd == 0 || _nr > nr_nd) {
    if(nr_nd != 0) {
      free_mms_tm_nd();
    }
    nr_nd = _nr;

    sigma = (double*)calloc((nr_nd), sizeof(double));
    zitam1 = (double*)calloc((nr_nd), sizeof(double));
    zita = (double*)calloc((nr_nd), sizeof(double));
    alphas = (double*)calloc((nr_nd), sizeof(double));
    betas = (double*)calloc((nr_nd), sizeof(double));

    ps_qmms = (spinor*)calloc(N*(2*nr_nd)+1,sizeof(spinor));
    ps_mms_solver = (spinor**)calloc((2*nr_nd)+1,sizeof(spinor*));

    for(int i = 0; i < 2*nr_nd; i++) {
      ps_mms_solver[i]=(spinor*)(((unsigned long int)(ps_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*N;
    }
    ini_mms_nd = 1;
  }
}

static void free_mms_tm_nd() {
  free(sigma);
  free(zitam1);
  free(zita);
  free(alphas);
  free(betas);
  free(ps_qmms);
  free(ps_mms_solver);
  nr_nd = 0;
  ini_mms_nd = 0;
  return;
}
