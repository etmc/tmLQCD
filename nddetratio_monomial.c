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
#include "linsolve.h"
#include "solver/solver.h"
#include "deriv_Sb.h"
#include "tm_operators.h"
#include "chebyshev_polynomial.h"
#include "Nondegenerate_Matrix.h"
#include "Hopping_Matrix.h"
#include "phmc.h"
#include "boundary.h"
#include "gamma.h"
#include "Nondegenerate_Matrix.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "reweighting_factor_nd.h"
#include "monomial.h"
#include "nddetratio_monomial.h"



double nddetratio_acc(const int id) {
  int iter;
  monomial * mnl = &monomial_list[id];

  
  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  boundary(mnl->kappa);

  iter = cg_her_nd(g_spinor_field[2], g_spinor_field[3], mnl->pf, mnl->pf2,
		   mnl->maxiter, mnl->accprec, g_relative_precision_flag, 
		   VOLUME/2, &Q_Qdagger_ND);
  QdaggerNon_degenerate(g_spinor_field[0], g_spinor_field[1],
			g_spinor_field[2], g_spinor_field[3]);

  g_mubar = mnl->mubar2;
  g_epsbar = mnl->epsbar2;
  boundary(mnl->kappa2);

  QNon_degenerate(g_spinor_field[2], g_spinor_field[3],
		  g_spinor_field[0], g_spinor_field[1]);
  
  mnl->energy1  = scalar_prod_r(mnl->pf , g_spinor_field[2], VOLUME/2, 1);
  mnl->energy1 += scalar_prod_r(mnl->pf2, g_spinor_field[3], VOLUME/2, 1);

  return(mnl->energy1 - mnl->energy0);
}
