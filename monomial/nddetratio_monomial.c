/***********************************************************************
 *
 * Copyright (C) 2008 Thomas Chiarappa, Carsten Urbach
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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/solver.h"
#include "deriv_Sb.h"
#include "operator/tm_operators.h"
#include "chebyshev_polynomial.h"
#include "operator/tm_operators_nd.h"
#include "operator/Hopping_Matrix.h"
#include "phmc.h"
#include "boundary.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "gamma.h"
#include "operator/tm_operators_nd.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "gettime.h"
#include "reweighting_factor_nd.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "nddetratio_monomial.h"
#include "DDalphaAMG_interface.h"
#include "operator/tm_operators_nd_32.h"

int init_nddetratio(rational_t * rat) {
  rat->np = 1;
  if((rat->mu = (double*)malloc(rat->np*sizeof(double))) == NULL)  {
    fprintf(stderr, "Could not allocate memory for coefficients in init_rational\n");
    return(-2);
  }
  rat->mu[0] = 0;
}

double nddetratio_acc(const int id, hamiltonian_field_t * const hf) {
  int iter;
  monomial * mnl = &monomial_list[id];
  
  double kappa_tmp;
  double csw_tmp;
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  matrix_mult_nd Q_pm_ndpsi = Qtm_pm_ndpsi, Q_dagger_ndpsi = Qtm_dagger_ndpsi, Q_ndpsi = Qtm_ndpsi;
  matrix_mult_nd32 Q_pm_ndpsi_32 = Qtm_pm_ndpsi_32;
  
  kappa_tmp=g_kappa;
  csw_tmp=g_c_sw;

  g_kappa=mnl->kappa;
  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  boundary(mnl->kappa);
  if(mnl->type == NDDETRATIO) {
    // NDDETRATIO does not have clover term so we need to set c_sw to zero otherwise QUDA will use it
    g_c_sw = 0;
  }
  if(mnl->type == NDCLOVERDETRATIO) {
    Q_pm_ndpsi = Qsw_pm_ndpsi;
    Q_pm_ndpsi_32 = Qsw_pm_ndpsi_32;
    Q_dagger_ndpsi = Qsw_dagger_ndpsi;
    Q_ndpsi = Qsw_ndpsi;
    init_sw_fields();
    sw_term((const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
  }
  if( mnl->solver == MG ) {
    iter = MG_solver_nd(mnl->w_fields[2], mnl->w_fields[3], mnl->pf, mnl->pf2,
                        mnl->accprec, mnl->maxiter, g_relative_precision_flag, 
                        VOLUME/2, g_gauge_field, Q_ndpsi);
  } else {
    // iter = cg_her_nd(mnl->w_fields[0], mnl->w_fields[1], mnl->pf, mnl->pf2,
    //                  mnl->maxiter, mnl->accprec, g_relative_precision_flag, 
    //                  VOLUME/2, Q_pm_ndpsi);
    mnl->solver_params.max_iter = mnl->maxiter;
    mnl->solver_params.squared_solver_prec = mnl->accprec;
    mnl->solver_params.rel_prec = g_relative_precision_flag;               
    mnl->solver_params.sdim = VOLUME/2;
    mnl->solver_params.no_shifts = mnl->rat.np; // set to 1 by init_nddetratio
    mnl->solver_params.shifts = mnl->rat.mu; // set to 0 by init_nddetratio
    mnl->solver_params.type = mnl->solver;
    mnl->solver_params.M_ndpsi = Q_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = Q_pm_ndpsi_32; 
    spinor ** out0 = (spinor**)malloc((mnl->rat.np)*sizeof(spinor*));
    spinor ** out1 = (spinor**)malloc((mnl->rat.np)*sizeof(spinor*));
    for (int j=0;j<mnl->rat.np;j++) {
      out0[j] = mnl->w_fields[0];
      out1[j] = mnl->w_fields[1];
    }
    iter = solve_mms_nd(out0, out1, mnl->pf, mnl->pf2,
                        &(mnl->solver_params));
    for (int j=0;j<mnl->rat.np;j++) {
      out0[j] = NULL;
      out1[j] = NULL;
    }
    free(out0); free(out1);
    tm_stopwatch_push(&g_timers, "Q_dagger_ndpsi", "");
    Q_dagger_ndpsi(mnl->w_fields[2], mnl->w_fields[3],
                   mnl->w_fields[0], mnl->w_fields[1]);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
  }

  g_kappa=mnl->kappa2;
  g_mubar = mnl->mubar2;
  g_epsbar = mnl->epsbar2;
  boundary(mnl->kappa2);

  if(mnl->type == NDCLOVERDETRATIO) {
    sw_term((const su3**) hf->gaugefield, mnl->kappa2, mnl->c_sw); 
    sw_invert_nd(mnl->mubar2*mnl->mubar2 - mnl->epsbar2*mnl->epsbar2);
  }
  tm_stopwatch_push(&g_timers, "Q_ndpsi", "");
  Q_ndpsi(mnl->w_fields[0], mnl->w_fields[1],
            mnl->w_fields[2], mnl->w_fields[3]);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
 
  tm_stopwatch_push(&g_timers, "energy1", ""); 
  mnl->energy1  = scalar_prod_r(mnl->pf , mnl->w_fields[0], VOLUME/2, 1);
  mnl->energy1 += scalar_prod_r(mnl->pf2, mnl->w_fields[1], VOLUME/2, 1);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called nddetratio_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy0 - mnl->energy1);
    }
  }
  
  g_kappa = kappa_tmp;
  g_c_sw = csw_tmp;
  
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return(mnl->energy1 - mnl->energy0);
}
