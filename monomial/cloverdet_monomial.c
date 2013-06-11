/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
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
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "ranlxd.h"
#include "sse.h"
#include "start.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "gamma.h"
#include "operator/tm_operators.h"
#include "operator/Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "solver/solver.h"
#include "operator/clover_leaf.h"
#include "read_input.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "monomial/monomial.h"
#include "operator/clovertm_operators.h"
#include "cloverdet_monomial.h"

/* think about chronological solver ! */

void cloverdet_derivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();
  for(int i = 0; i < VOLUME; i++) { 
    for(int mu = 0; mu < 4; mu++) { 
      _su3_zero(swm[i][mu]);
      _su3_zero(swp[i][mu]);
    }
  }

  mnl->forcefactor = 1.;
  /*********************************************************************
   * 
   * even/odd version 
   *
   * This a term is det(\hat Q^2(\mu))
   *
   *********************************************************************/
  
  g_mu = mnl->mu;
  g_mu3 = mnl->rho;
  boundary(mnl->kappa);
  
  // we compute the clover term (1 + T_ee(oo)) for all sites x
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  // we invert it for the even sites only
  sw_invert(EE, mnl->mu);
  
  if(mnl->solver != CG && g_proc_id == 0) {
    fprintf(stderr, "Bicgstab currently not implemented, using CG instead! (cloverdet_monomial.c)\n");
  }
  
  // Invert Q_{+} Q_{-}
  // X_o -> w_fields[1]
  chrono_guess(mnl->w_fields[1], mnl->pf, mnl->csg_field, mnl->csg_index_array,
	       mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
  mnl->iter1 += cg_her(mnl->w_fields[1], mnl->pf, mnl->maxiter, mnl->forceprec, 
		       g_relative_precision_flag, VOLUME/2, mnl->Qsq);
  chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);
  
  // Y_o -> w_fields[0]
  mnl->Qm(mnl->w_fields[0], mnl->w_fields[1]);
  
  // apply Hopping Matrix M_{eo}
  // to get the even sites of X_e
  H_eo_sw_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1, mnl->mu);
  // \delta Q sandwitched by Y_o^\dagger and X_e
  deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor); 
  
  // to get the even sites of Y_e
  H_eo_sw_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1, mnl->mu);
  // \delta Q sandwitched by Y_e^\dagger and X_o
  // uses the gauge field in hf and changes the derivative fields in hf
  deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf, mnl->forcefactor);
  
  // here comes the clover term...
  // computes the insertion matrices for S_eff
  // result is written to swp and swm
  // even/even sites sandwiched by gamma_5 Y_e and gamma_5 X_e
  sw_spinor(EE, mnl->w_fields[2], mnl->w_fields[3], mnl->forcefactor);
  
  // odd/odd sites sandwiched by gamma_5 Y_o and gamma_5 X_o
  sw_spinor(OO, mnl->w_fields[0], mnl->w_fields[1], mnl->forcefactor);
  
  // compute the contribution for the det-part
  // we again compute only the insertion matrices for S_det
  // the result is added to swp and swm
  // even sites only!
  sw_deriv(EE, mnl->mu);
  
  // now we compute
  // finally, using the insertion matrices stored in swm and swp
  // we compute the terms F^{det} and F^{sw} at once
  // uses the gaugefields in hf and changes the derivative field in hf
  sw_all(hf, mnl->kappa, mnl->c_sw);

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}


void cloverdet_heatbath(const int id, hamiltonian_field_t * const hf) {

  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();

  g_mu = mnl->mu;
  g_mu3 = mnl->rho;
  g_c_sw = mnl->c_sw;
  boundary(mnl->kappa);
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;

  init_sw_fields();
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EE, mnl->mu);

  random_spinor_field_eo(mnl->w_fields[0], mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME/2, 1);
  
  mnl->Qp(mnl->pf, mnl->w_fields[0]);
  chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called cloverdet_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}


double cloverdet_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  int save_sloppy = g_sloppy_precision_flag;
  double atime, etime;
  atime = gettime();

  g_mu = mnl->mu;
  g_mu3 = mnl->rho;
  g_c_sw = mnl->c_sw;
  boundary(mnl->kappa);

  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EE, mnl->mu);

  chrono_guess(mnl->w_fields[0], mnl->pf, mnl->csg_field, mnl->csg_index_array,
	       mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
  g_sloppy_precision_flag = 0;
  mnl->iter0 = cg_her(mnl->w_fields[0], mnl->pf, mnl->maxiter, mnl->accprec,  
		      g_relative_precision_flag, VOLUME/2, mnl->Qsq); 
  mnl->Qm(mnl->w_fields[0], mnl->w_fields[0]);
  
  g_sloppy_precision_flag = save_sloppy;
  /* Compute the energy contr. from first field */
  mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME/2, 1);

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called cloverdet_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy1 - mnl->energy0);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}
