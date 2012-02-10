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
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "gamma.h"
#include "tm_operators.h"
#include "Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "solver/solver.h"
#include "clover_leaf.h"
#include "read_input.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "monomial.h"
#include "clover.h"
#include "cloverdet_monomial.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

/* think about chronological solver ! */

void cloverdet_derivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];

  /* This factor 2 a missing factor 2 in trace_lambda */
  (*mnl).forcefactor = 2.;

  /*********************************************************************
   * 
   * even/odd version 
   *
   * This a term is det(\hat Q^2(\mu))
   *
   *********************************************************************/
  
  g_mu = 0.;
  boundary(mnl->kappa);
  
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EO);
  
  if(mnl->solver != CG && g_proc_id == 0) {
    fprintf(stderr, "Bicgstab currently not implemented, using CG instead! (cloverdet_monomial.c)\n");
  }
  
  /* Invert Q_{+} Q_{-} */
  /* X_o -> DUM_DERI+1 */
  chrono_guess(g_spinor_field[DUM_DERI+1], mnl->pf, mnl->csg_field, mnl->csg_index_array,
	       mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
  mnl->iter1 += cg_her(g_spinor_field[DUM_DERI+1], mnl->pf, mnl->maxiter, mnl->forceprec, 
		       g_relative_precision_flag, VOLUME/2, mnl->Qsq);
  chrono_add_solution(g_spinor_field[DUM_DERI+1], mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);
  
  /* Y_o -> DUM_DERI  */
  mnl->Qm(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  
  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X_e */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO);
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], hf); 
  
  /* to get the even sites of Y_e */
  H_eo_sw_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO);
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1], hf);
  
  /* here comes the clover term... */
  gamma5(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+2], VOLUME/2);
  sw_spinor(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
  
  /* compute the contribution for the det-part */
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  sw_spinor(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  
  sw_deriv(EO);
  sw_all(hf, mnl->kappa, mnl->c_sw);

  g_mu = g_mu1;
  boundary(g_kappa);

  return;
}


void cloverdet_heatbath(const int id, hamiltonian_field_t * const hf) {

  monomial * mnl = &monomial_list[id];
  g_mu = 0.;
  boundary(mnl->kappa);
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;

  init_sw_fields();
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EO);

  random_spinor_field(g_spinor_field[2], VOLUME/2, mnl->rngrepro);
  mnl->energy0 = square_norm(g_spinor_field[2], VOLUME/2, 1);
  
  mnl->Qp(mnl->pf, g_spinor_field[2]);
  chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
		      mnl->csg_N, &mnl->csg_n, VOLUME/2);

  g_mu = g_mu1;
  boundary(g_kappa);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverdet_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }
  return;
}


double cloverdet_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  int save_sloppy = g_sloppy_precision_flag;

  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert(EO);

  g_mu = 0.;
  boundary(mnl->kappa);
  
  chrono_guess(g_spinor_field[2], mnl->pf, mnl->csg_field, mnl->csg_index_array,
	       mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
  g_sloppy_precision_flag = 0;
  mnl->iter0 = cg_her(g_spinor_field[2], mnl->pf, mnl->maxiter, mnl->accprec,  
		      g_relative_precision_flag, VOLUME/2, mnl->Qsq); 
  mnl->Qm(g_spinor_field[2], g_spinor_field[2]);
  
  g_sloppy_precision_flag = save_sloppy;
  /* Compute the energy contr. from first field */
  mnl->energy1 = square_norm(g_spinor_field[2], VOLUME/2, 1);

  g_mu = g_mu1;
  boundary(g_kappa);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverdet_acc for id %d %d dH = %1.4e\n", 
	   id, mnl->even_odd_flag, mnl->energy1 - mnl->energy0);
  }
  return(mnl->energy1 - mnl->energy0);
}
