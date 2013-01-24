/***********************************************************************
 *
 * Copyright (C) 2013 Carsten Urbach
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
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/solver.h"
#include "deriv_Sb.h"
#include "init/init_chi_spinor_field.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_nd.h"
#include "operator/Hopping_Matrix.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "rational/rational.h"
#include "phmc.h"
#include "rat_monomial.h"

// DEBUG: we should use here the not normalised operators
//        and, therefore, rescale the coefficients in the rational
//        approximation!
//
// DEBUG: light CGMMS solver needs to be re-written for 
//        this here
//
// This monomial can only work for mnl->mu = 0!
// Hence, we have to ensure it in init_monomial

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void rat_derivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  solver_pm_t solver_pm;
  double atime, etime;
  atime = gettime();
  g_mu = 0;
  g_mu3 = 0.;
  boundary(mnl->kappa);

  if(mnl->type == CLOVERRAT) {
    for(int i = 0; i < VOLUME; i++) { 
      for(int mu = 0; mu < 4; mu++) { 
	_su3_zero(swm[i][mu]);
	_su3_zero(swp[i][mu]);
      }
    }
  
    // we compute the clover term (1 + T_ee(oo)) for all sites x
    sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    // we invert it for the even sites only
    sw_invert(EE, mnl->mu);
  }
  //mnl->forcefactor = mnl->EVMaxInv;
  mnl->forcefactor = 1.;

  solver_pm.max_iter = mnl->maxiter;
  solver_pm.squared_solver_prec = mnl->forceprec;
  solver_pm.no_shifts = mnl->rat.np;
  solver_pm.shifts = mnl->rat.mu;
  solver_pm.rel_prec = g_relative_precision_flag;
  solver_pm.type = CGMMS;
  solver_pm.M_psi = mnl->Qsq;
  solver_pm.sdim = VOLUME/2;
  // this generates all X_j,o (odd sites only) -> g_chi_up|dn_spinor_field
  //mnl->iter1 += cg_mms_tm(g_chi_up_spinor_field, mnl->pf,
  //		     &solver_pm);
  
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    mnl->Qp(mnl->w_fields[0], g_chi_up_spinor_field[j]);
    if(mnl->type == CLOVERRAT) {
      // apply Hopping Matrix M_{eo}
      // to get the even sites of X_e
      H_eo_sw_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[j], EO, -mnl->mu);
      // \delta Q sandwitched by Y_o^\dagger and X_e
      deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor); 
      
      // to get the even sites of Y_e
      H_eo_sw_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, mnl->mu);
      // \delta Q sandwitched by Y_e^\dagger and X_o
      // uses the gauge field in hf and changes the derivative fields in hf
      deriv_Sb(EO, mnl->w_fields[3], g_chi_up_spinor_field[j], hf, mnl->forcefactor);

    }
    else {
      /* apply Hopping Matrix M_{eo} */
      /* to get the even sites of X_e */
      H_eo_tm_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[j], EO, -1.);
      /* \delta Q sandwitched by Y_o^\dagger and X_e */
      deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor); 
      
      /* to get the even sites of Y_e */
      H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
      /* \delta Q sandwitched by Y_e^\dagger and X_o */
      deriv_Sb(EO, mnl->w_fields[3], g_chi_up_spinor_field[j], hf, mnl->forcefactor);
    }

    if(mnl->type == CLOVERRAT) {
      // even/even sites sandwiched by gamma_5 Y_e and gamma_5 X_e
      sw_spinor(EE, mnl->w_fields[2], mnl->w_fields[3], mnl->forcefactor);
  
      // odd/odd sites sandwiched by gamma_5 Y_o and gamma_5 X_o
      sw_spinor(OO, mnl->w_fields[0], g_chi_up_spinor_field[j], mnl->forcefactor);
    }
  }
  // trlog part does not depend on the normalisation
  // here we need another mechanism!
  if(mnl->type == CLOVERRAT) {
    sw_deriv(EE, mnl->mu);
    sw_all(hf, mnl->kappa, mnl->c_sw);
  }
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}


void rat_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  solver_pm_t solver_pm;
  double atime, etime;
  atime = gettime();
  // only for non-twisted operators
  g_mu = 0.;
  g_mu3 = 0.;
  mnl->mu = 0.;
  boundary(mnl->kappa);

  mnl->iter1 = 0;
  g_mu3 = 0.;
  if(mnl->type == CLOVERRAT) {
    init_sw_fields();
    sw_term((const su3**)hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert(EE, mnl->mu);
  }
  // we measure before the trajectory!
  if((mnl->rec_ev != 0) && (hf->traj_counter%mnl->rec_ev == 0)) {
    //if(mnl->type != NDCLOVERRAT) phmc_compute_ev(hf->traj_counter-1, id, &Qtm_pm_ndbipsi);
    //else phmc_compute_ev(hf->traj_counter-1, id, &Qsw_pm_ndbipsi);
  }

  // the Gaussian distributed random fields
  mnl->energy0 = 0.;
  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(mnl->pf, VOLUME/2, 1);

  // set solver parameters
  solver_pm.max_iter = mnl->maxiter;
  solver_pm.squared_solver_prec = mnl->accprec;
  solver_pm.no_shifts = mnl->rat.np;
  solver_pm.shifts = mnl->rat.nu;
  solver_pm.type = CGMMS;
  solver_pm.M_psi = mnl->Qsq;
  solver_pm.sdim = VOLUME/2;
  solver_pm.rel_prec = g_relative_precision_flag;
  //mnl->iter0 = cg_mms_tm_nd(g_chi_up_spinor_field, mnl->pf,
  //			     &solver_pm);

  assign(mnl->w_fields[2], mnl->pf, VOLUME/2);

  // apply C to the random field to generate pseudo-fermion fields
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    // Q^m - i nu_j
    mnl->Qp(g_chi_up_spinor_field[mnl->rat.np], g_chi_up_spinor_field[j]);
    assign_add_mul(g_chi_up_spinor_field[mnl->rat.np], g_chi_up_spinor_field[j], -I*mnl->rat.nu[j], VOLUME/2);
    assign_add_mul(mnl->pf, g_chi_up_spinor_field[mnl->rat.np], I*mnl->rat.rnu[j], VOLUME/2);
  }

  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) { 
      printf("called ndrat_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}


double rat_acc(const int id, hamiltonian_field_t * const hf) {
  solver_pm_t solver_pm;
  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();
  // only for non-twisted operators
  g_mu = 0.;
  g_mu3 = 0.;
  boundary(mnl->kappa);
  if(mnl->type == CLOVERRAT) {
    sw_term((const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert(EE, mnl->mu);
  }
  mnl->energy1 = 0.;

  solver_pm.max_iter = mnl->maxiter;
  solver_pm.squared_solver_prec = mnl->accprec;
  solver_pm.no_shifts = mnl->rat.np;
  solver_pm.shifts = mnl->rat.mu;
  solver_pm.type = CGMMS;
  solver_pm.M_psi = mnl->Qsq;
  solver_pm.sdim = VOLUME/2;
  solver_pm.rel_prec = g_relative_precision_flag;
  //mnl->iter0 += cg_mms_tm_nd(g_chi_up_spinor_field, mnl->pf,
  //		     &solver_pm);

  // apply R to the pseudo-fermion fields
  assign(mnl->w_fields[0], mnl->pf, VOLUME/2);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(mnl->w_fields[0], g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }

  mnl->energy1 = scalar_prod_r(mnl->pf, mnl->w_fields[0], VOLUME/2, 1);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 0) { // shoud be 3
      printf("called ndrat_acc for id %d dH = %1.10e\n", id, mnl->energy1 - mnl->energy0);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}


