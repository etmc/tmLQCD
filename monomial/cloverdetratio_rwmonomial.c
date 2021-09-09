/***********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "gamma.h"
#include "operator/tm_operators.h"
#include "operator/Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "solver/solver.h"
#include "solver/monomial_solve.h"
#include "read_input.h"
#include "operator/clovertm_operators.h"
#include "operator/clovertm_operators_32.h"
#include "operator/clover_leaf.h"
#include "monomial/monomial.h"
#include "boundary.h"
#include "cloverdetratio_rwmonomial.h"
#include "expo.h"
#include "xchange/xchange.h"
#include "init/init_gauge_tmp.h"
#include "DDalphaAMG_interface.h"

double cloverdetratio_rwacc(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  int save_sloppy = g_sloppy_precision_flag;

  g_mu = mnl->mu2;
  boundary(mnl->kappa2);

  init_sw_fields();
  sw_term( (const su3**) hf->gaugefield, mnl->kappa2, mnl->c_sw); 
  sw_invert(EE, mnl->mu2);
  g_mu3 = 0.;
  mnl->Qp(mnl->w_fields[1], mnl->pf);

  g_mu3 = 0.;
  g_mu = mnl->mu;
  boundary(mnl->kappa);
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EE, mnl->mu);

  chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array, 
	       mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_plus_psi);
  g_sloppy_precision_flag = 0;    
  if( mnl->solver == MG || mnl->solver == BICGSTAB ) {
    mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params, 
                                   mnl->maxiter, mnl->accprec,
				   g_relative_precision_flag, VOLUME/2, mnl->Qp, mnl->solver);
  } else {
    mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params, mnl->maxiter, mnl->accprec,
				   g_relative_precision_flag, VOLUME/2, mnl->Qsq, mnl->solver);
    tm_stopwatch_push(&g_timers);
    mnl->Qm(mnl->w_fields[0], mnl->w_fields[0]);
    tm_stopwatch_pop(&g_timers, 0, 1, "", "Qm");
  }

  g_sloppy_precision_flag = save_sloppy;

  /* Compute the energy contr. from second field */
  tm_stopwatch_push(&g_timers);
  mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME/2, 1);
  tm_stopwatch_pop(&g_timers, 0, 1, "", "energy1");

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called cloverdetratio_rwacc for id %d dH = %1.10e\n", 
	     id, mnl->energy1 - mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return(mnl->energy1 - mnl->energy0);
}
