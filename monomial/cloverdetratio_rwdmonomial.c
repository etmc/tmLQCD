/***********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
 *
 * Adaption with a more efficient inversion by Georg Bergner 2016
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
# include<config.h>
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
#include "cloverdetratio_rwdmonomial.h"


double cloverdetratio_rwaccd(const int id, hamiltonian_field_t * const hf) {
	monomial * mnl = &monomial_list[id];
	  int save_sloppy = g_sloppy_precision_flag;
	  double atime, etime;
	  atime = gettime();


	  g_mu = mnl->mu2;
	  boundary(mnl->kappa2);

	  init_sw_fields();
	  sw_term( (const su3**) hf->gaugefield, mnl->kappa2, mnl->c_sw);
	  sw_invert(EE, mnl->mu2);
	  g_mu3 = 0.;

	  mnl->Qsq(mnl->w_fields[1], mnl->pf);
	  assign(mnl->w_fields[0],mnl->pf,VOLUME/2);

	  g_mu3 = 0.;
	  g_mu = mnl->mu;
	  boundary(mnl->kappa);
	  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw);
	  sw_invert(EE, mnl->mu);

	  /*chrono_guess(mnl->pf, mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
		       mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_plus_psi);*/
	  g_sloppy_precision_flag = 0;
	  mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params, mnl->maxiter, mnl->accprec,
			      g_relative_precision_flag, VOLUME/2, mnl->Qsq, mnl->solver);


	  g_sloppy_precision_flag = save_sloppy;

  /* Compute the energy contr. from second field */
  mnl->energy1 = scalar_prod_r(mnl->pf,mnl->w_fields[0], VOLUME/2, 1);

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial rwacc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called cloverdetratio_rwacc for id %d dH = %1.10e\n", 
	     id, mnl->energy1 - mnl->energy0);
  	  printf("Parameters: kappa=%e, kappa2=%e; mu=%e, mu2=%e, vprod=%e\n",mnl->kappa,mnl->kappa2,mnl->mu,mnl->mu2,mnl->energy1);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}
