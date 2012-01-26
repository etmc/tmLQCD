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
#include <io/spinor.h>
#include "gamma.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "cg_mms_tm.h"
#include "solver_field.h"
#include <io/params.h>

static spinor * xs_qmms;
static spinor * ps_qmms;
static spinor ** xs_mms_solver;
static spinor ** ps_mms_solver;
static double * sigma;
static double * zitam1, * zita;
static double * alphas, * betas;

extern int index_start;

void init_mms_tm(const int nr);


/* P output = solution , Q input = source */
int cg_mms_tm(spinor * const P, spinor * const Q, const int max_iter, 
	      double eps_sq, const int rel_prec, const int N, matrix_mult f) {

  static double normsq, pro, err, alpha_cg = 1., beta_cg = 0., squarenorm;
  int iteration, im, append = 0;
  char filename[100];
  static double gamma, alpham1;
  int const cg_mms_default_precision = 32;
  double tmp_mu = g_mu;
  WRITER * writer = NULL;
  paramsInverterInfo *inverterInfo = NULL;
  paramsPropagatorFormat *propagatorFormat = NULL;
  spinor * temp_save; //used to save all the masses
  spinor ** solver_field = NULL;
  const int nr_sf = 5;

  init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);
  init_mms_tm(g_no_extra_masses);

  /* currently only implemented for P=0 */
  zero_spinor_field(P, N);
  /*  Value of the bare MMS-masses (\mu^2 - \mu_0^2) */
  for(im = 0; im < g_no_extra_masses; im++) {
    sigma[im] = g_extra_masses[im]*g_extra_masses[im] - g_mu*g_mu;
    assign(xs_mms_solver[im], P, N);
    assign(ps_mms_solver[im], Q, N);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }

  squarenorm = square_norm(Q, N, 1);
  assign(solver_field[0], P, N);
/*   normsp = square_norm(P, N, 1); */

  /* initialize residue r and search vector p */
/*   if(normsp == 0){ */
  /* currently only implemented for P=0 */
  if(1) {
    /* if a starting solution vector equal to zero is chosen */
    assign(solver_field[1], Q, N);
    assign(solver_field[2], Q, N);
    normsq = square_norm(Q, N, 1);
  }
  else{
    /* if a starting solution vector different from zero is chosen */
    f(solver_field[3], solver_field[0]);

    diff(solver_field[1], Q, solver_field[3], N);
    assign(solver_field[2], solver_field[1], N);
    normsq = square_norm(solver_field[2], N, 1);
  }

  /* main loop */
  for(iteration = 0; iteration < max_iter; iteration++) {

    /*   Q^2*p and then (p,Q^2*p)  */
    f(solver_field[4], solver_field[2]);
    pro = scalar_prod_r(solver_field[2], solver_field[4], N, 1);

    /* For the update of the coeff. of the shifted pol. we need alpha_cg(i-1) and alpha_cg(i).
       This is the reason why we need this double definition of alpha */
    alpham1 = alpha_cg;

    /* Compute alpha_cg(i+1) */
    alpha_cg = normsq/pro;
    for(im = 0; im < g_no_extra_masses; im++) {

      /* Now gamma is a temp variable that corresponds to zita(i+1) */ 
      gamma = zita[im]*alpham1/(alpha_cg*beta_cg*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alpha_cg));

      /* Now zita(i-1) is put equal to the old zita(i) */
      zitam1[im] = zita[im];
      /* Now zita(i+1) is updated */
      zita[im] = gamma;
      /* Update of alphas(i) = alpha_cg(i)*zita(i+1)/zita(i) */ 
      alphas[im] = alpha_cg*zita[im]/zitam1[im];
      /* Compute xs(i+1) = xs(i) + alphas(i)*ps(i) */
      assign_add_mul_r(xs_mms_solver[im], ps_mms_solver[im], alphas[im], N); 
    }

    /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
    assign_add_mul_r(solver_field[0], solver_field[2],  alpha_cg, N);
    /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Qp_i   */
    assign_add_mul_r(solver_field[1], solver_field[4], -alpha_cg, N);

    /* Check whether the precision eps_sq is reached */

    err = square_norm(solver_field[1], N, 1);
    if(g_debug_level > 2 && g_proc_id == g_stdio_proc) {
      printf("CGMMS iteration: %d residue: %g\n", iteration, err); fflush( stdout );
    }

    if( ((err <= eps_sq) && (rel_prec == 0)) ||
      ((err <= eps_sq*squarenorm) && (rel_prec == 1)) ) {

      assign(P, solver_field[0], N);
      f(solver_field[2], P);
      diff(solver_field[3], solver_field[2], Q, N);
      err = square_norm(solver_field[3], N, 1);
      if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
        printf("# CG MMS true residue at final iteration (%d) was %g.\n", iteration, err); 
        fflush( stdout);
      }
      g_sloppy_precision = 0;
      g_mu = tmp_mu;

      /* save all the results of (Q^dagger Q)^(-1) \gamma_5 \phi */
      /* here ... */
      /* when im == -1 save the base mass*/
      for(im = -1; im < g_no_extra_masses; im++) {
        if(im==-1) {
          temp_save=solver_field[0];
        } else {
          temp_save=xs_mms_solver[im];
        }

        if(SourceInfo.type != 1) {
          if (PropInfo.splitted) {
            sprintf(filename, "%s.%.4d.%.2d.%.2d.cgmms.%.2d.inverted", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, SourceInfo.ix, im+1);
          } else {
            sprintf(filename, "%s.%.4d.%.2d.cgmms.%.2d.inverted", SourceInfo.basename, SourceInfo.nstore, SourceInfo.t, im+1);
          }
        }
        else {
          sprintf(filename, "%s.%.4d.%.5d.cgmms.%.2d.0", SourceInfo.basename, SourceInfo.nstore, SourceInfo.sample, im+1);
        }
        if(g_kappa != 0) {
          mul_r(temp_save, (2*g_kappa)*(2*g_kappa), temp_save, N);
        }

        append = !PropInfo.splitted;

        construct_writer(&writer, filename, append);

        if (PropInfo.splitted || SourceInfo.ix == index_start) {
          //Create the inverter info NOTE: always set to TWILSON=12 and 1 flavour (to be adjusted)
          inverterInfo = construct_paramsInverterInfo(err, iteration+1, 12, 1);
          if (im == -1) {
            inverterInfo->cgmms_mass = inverterInfo->mu;
          } else {
            inverterInfo->cgmms_mass = g_extra_masses[im]/(2 * inverterInfo->kappa);
          }
          write_spinor_info(writer, PropInfo.format, inverterInfo, append);
          //Create the propagatorFormat NOTE: always set to 1 flavour (to be adjusted)
          propagatorFormat = construct_paramsPropagatorFormat(cg_mms_default_precision, 1);
          write_propagator_format(writer, propagatorFormat);
          free(inverterInfo);
          free(propagatorFormat);
        }
        convert_lexic_to_eo(solver_field[2], solver_field[1], temp_save);
        write_spinor(writer, &solver_field[2], &solver_field[1], 1, cg_mms_default_precision);
        destruct_writer(writer);
      }
      finalize_solver(solver_field, nr_sf);
      return(iteration+1);
    }

    /* Compute beta_cg(i+1) = (r(i+1),r(i+1))/(r(i),r(i))
       Compute p(i+1) = r(i+1) + beta(i+1)*p(i)  */
    beta_cg = err/normsq;
    assign_mul_add_r(solver_field[2], beta_cg, solver_field[1], N);
    normsq = err;

    /* Compute betas(i+1) = beta_cg(i)*(zita(i+1)*alphas(i))/(zita(i)*alpha_cg(i))
       Compute ps(i+1) = zita(i+1)*r(i+1) + betas(i+1)*ps(i)  */
    for(im = 0; im < g_no_extra_masses; im++) {
      betas[im] = beta_cg*zita[im]*alphas[im]/(zitam1[im]*alpha_cg);
      assign_mul_add_mul_r(ps_mms_solver[im], solver_field[1], betas[im], zita[im], N);
    }
  }
  assign(P, solver_field[0], N);
  g_sloppy_precision = 0;
  finalize_solver(solver_field, nr_sf);
  return(-1);
}


void init_mms_tm(const int nr) {
  static int ini = 0;
  int i;
  if(ini == 0) {

    sigma = (double*)calloc((nr), sizeof(double));
    zitam1 = (double*)calloc((nr), sizeof(double));
    zita = (double*)calloc((nr), sizeof(double));
    alphas = (double*)calloc((nr), sizeof(double));
    betas = (double*)calloc((nr), sizeof(double));

#if (defined SSE2 || defined SSE)
    xs_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr)+1,sizeof(spinor));
    xs_mms_solver = (spinor**)calloc((nr)+1,sizeof(spinor*));
    
    ps_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr)+1,sizeof(spinor));
    ps_mms_solver = (spinor**)calloc((nr)+1,sizeof(spinor*));

    for(i = 0; i < nr; i++) {
      xs_mms_solver[i]=(spinor*)(((unsigned long int)(xs_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*VOLUMEPLUSRAND;
      ps_mms_solver[i]=(spinor*)(((unsigned long int)(ps_qmms)+ALIGN_BASE)&~ALIGN_BASE) + i*VOLUMEPLUSRAND;
    }
#else
    xs_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr),sizeof(spinor));
    xs_mms_solver = (spinor**)calloc((nr),sizeof(spinor*));

    ps_qmms = (spinor*)calloc(VOLUMEPLUSRAND*(nr),sizeof(spinor));
    ps_mms_solver = (spinor**)calloc((nr),sizeof(spinor*));

    for(i = 0; i < nr; i++) {
      xs_mms_solver[i] = xs_qmms + i*VOLUMEPLUSRAND;
      ps_mms_solver[i] = ps_qmms + i*VOLUMEPLUSRAND;
    }
#endif
    ini=1;
  }
}
