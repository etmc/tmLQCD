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
# include<config.h>
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
#include "gamma.h"
#include "operator/tm_operators_nd.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "gettime.h"
#include "reweighting_factor_nd.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "nddetratio_monomial.h"



double nddetratio_acc(const int id, hamiltonian_field_t * const hf) {
  int iter;
  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();
  
  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  boundary(mnl->kappa);

  iter = cg_her_nd(mnl->w_fields[0], mnl->w_fields[1], mnl->pf, mnl->pf2,
		   mnl->maxiter, mnl->accprec, g_relative_precision_flag, 
		   VOLUME/2, &Qtm_pm_ndpsi);
  Qtm_dagger_ndpsi(mnl->w_fields[2], mnl->w_fields[3],
			mnl->w_fields[0], mnl->w_fields[1]);

  g_mubar = mnl->mubar2;
  g_epsbar = mnl->epsbar2;
  boundary(mnl->kappa2);

  Qtm_ndpsi(mnl->w_fields[0], mnl->w_fields[1],
		  mnl->w_fields[2], mnl->w_fields[3]);
  
  mnl->energy1  = scalar_prod_r(mnl->pf , mnl->w_fields[0], VOLUME/2, 1);
  mnl->energy1 += scalar_prod_r(mnl->pf2, mnl->w_fields[1], VOLUME/2, 1);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called nddetratio_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy0 - mnl->energy1);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}
