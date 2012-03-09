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
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "ranlxd.h"
#include "start.h"
#include "linalg_eo.h"
#include "linsolve.h"
#include "deriv_Sb.h"
#include "deriv_Sb_D_psi.h"
#include "gamma.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "solver/solver.h"
#include "read_input.h"
#include "smearing/stout.h"
#include "clover.h"
#include "clover_leaf.h"

#include "monomial.h"
#include "boundary.h"
#include "cloverdetratio_monomial.h"

/* think about chronological solver ! */

void cloverdetratio_derivative(const int no, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[no];

  /* This factor 2* a missing factor 2 in trace_lambda */
  mnl->forcefactor = 2.;

  /*********************************************************************
   *
   *  this is being run in case there is even/odd preconditioning
   * 
   * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
   * mu1 and mu2 are set according to the monomial
   *
   *********************************************************************/
  /* First term coming from the second field */
  /* Multiply with W_+ */
  g_mu = mnl->mu;
  g_mu3 = mnl->rho2; //rho2
  boundary(mnl->kappa);

  // we compute the clover term (1 + T_ee(oo)) for all sites x
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  // we invert it for the even sites only including mu
  sw_invert(EE, mnl->mu);
  
  if(mnl->solver != CG) {
    fprintf(stderr, "Bicgstab currently not implemented, using CG instead! (detratio_monomial.c)\n");
  }
  
  mnl->Qp(g_spinor_field[DUM_DERI+2], mnl->pf);
  g_mu3 = mnl->rho; // rho1

  /* Invert Q_{+} Q_{-} */
  /* X_W -> DUM_DERI+1 */
  chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], mnl->csg_field, 
	       mnl->csg_index_array, mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
  mnl->iter1 += cg_her(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], mnl->maxiter, 
		       mnl->forceprec, g_relative_precision_flag, VOLUME/2, mnl->Qsq);
  chrono_add_solution(g_spinor_field[DUM_DERI+1], mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);
  /* Y_W -> DUM_DERI  */
  mnl->Qm(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  
  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EE, -mnl->mu);
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], hf); 
  
  /* to get the even sites of Y */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EE, mnl->mu);
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1], hf); 

  // here comes the clover term...
  // computes the insertion matrices for S_eff
  // result is written to swp and swm
  // even/even sites sandwiched by gamma_5 Y_e and gamma_5 X_e  
  gamma5(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+2], VOLUME/2);
  sw_spinor(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
  
  // odd/odd sites sandwiched by gamma_5 Y_o and gamma_5 X_o
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  sw_spinor(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);

  g_mu3 = mnl->rho2; // rho2
  
  /* Second term coming from the second field */
  /* The sign is opposite!! */
  mul_r(g_spinor_field[DUM_DERI], -1., mnl->pf, VOLUME/2);
  
  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EE, -mnl->mu);
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], hf); 
  
  /* to get the even sites of Y */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EE, mnl->mu);
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1], hf);

  // here comes the clover term...
  // computes the insertion matrices for S_eff
  // result is written to swp and swm
  // even/even sites sandwiched by gamma_5 Y_e and gamma_5 X_e
  gamma5(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+2], VOLUME/2);
  sw_spinor(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
  
  // odd/odd sites sandwiched by gamma_5 Y_o and gamma_5 X_o
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  sw_spinor(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);

  sw_all(hf, mnl->kappa, mnl->c_sw);
  
  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);

  return;
}


void cloverdetratio_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];

  g_mu = mnl->mu;
  g_c_sw = mnl->c_sw;
  boundary(mnl->kappa);
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;
  
  init_sw_fields();
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EE, mnl->mu);

  random_spinor_field(g_spinor_field[4], VOLUME/2, mnl->rngrepro);
  mnl->energy0  = square_norm(g_spinor_field[4], VOLUME/2, 1);
  
  g_mu3 = mnl->rho;
  mnl->Qp(g_spinor_field[3], g_spinor_field[4]);
  g_mu3 = mnl->rho2;
  zero_spinor_field(mnl->pf,VOLUME/2);

  mnl->iter0 = cg_her(mnl->pf, g_spinor_field[3], mnl->maxiter, mnl->accprec,  
		      g_relative_precision_flag, VOLUME/2, mnl->Qsq); 

  chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);
  mnl->Qm(mnl->pf, mnl->pf);

  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverdetratio_heatbath for id %d %d energy %f\n", id, mnl->even_odd_flag, mnl->energy0);
  }
  g_mu3 = 0.;
  g_mu = g_mu1;
  boundary(g_kappa);
  return;
}

double cloverdetratio_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  int save_sloppy = g_sloppy_precision_flag;

  g_mu = mnl->mu;
  boundary(mnl->kappa);
  
  g_mu3 = mnl->rho2;
  mnl->Qp(g_spinor_field[DUM_DERI+5], mnl->pf);
  g_mu3 = mnl->rho;

  chrono_guess(g_spinor_field[3], g_spinor_field[DUM_DERI+5], mnl->csg_field, mnl->csg_index_array, 
	       mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_plus_psi);
  g_sloppy_precision_flag = 0;    
  mnl->iter0 += cg_her(g_spinor_field[3], g_spinor_field[DUM_DERI+5], mnl->maxiter, mnl->accprec,  
		      g_relative_precision_flag, VOLUME/2, mnl->Qsq);
  mnl->Qm(g_spinor_field[3], g_spinor_field[3]);

  g_sloppy_precision_flag = save_sloppy;

  /* Compute the energy contr. from second field */
  mnl->energy1 = square_norm(g_spinor_field[3], VOLUME/2, 1);

  g_mu = g_mu1;
  g_mu3 = 0.;
  boundary(g_kappa);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverdetratio_acc for id %d %d dH = %1.4e\n", 
	   id, mnl->even_odd_flag, mnl->energy1 - mnl->energy0);
  }
  return(mnl->energy1 - mnl->energy0);
}
