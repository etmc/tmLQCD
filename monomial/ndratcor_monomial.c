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
#include "ndrat_monomial.h"
#include "ndratcor_monomial.h"

// computes ||(1 - C^dagger R C) phi||
void check_C_ndpsi(spinor * const k_up, spinor * const k_dn,
		   spinor * const l_up, spinor * const l_dn,
		   const int id, hamiltonian_field_t * const hf,
		   solver_pm_t * solver_pm);

// applies (Q^2 R^2 -1) phi
double apply_Z_ndpsi(spinor * const k_up, spinor * const k_dn,
		     spinor * const l_up, spinor * const l_dn,
		     const int id, hamiltonian_field_t * const hf,
		     solver_pm_t * solver_pm);



void ndratcor_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  solver_pm_t solver_pm;
  double atime, etime, delta;
  spinor * up0, * dn0, * up1, * dn1, * tup, * tdn;
  double coefs[6] = {1./4., -3./32., 7./122., -77./2048., 231./8192., -1463./65536.};
  atime = gettime();
  nd_set_global_parameter(mnl);
  g_mu3 = 0.;
  mnl->iter0 = 0;
  if(mnl->type == NDCLOVERRATCOR) {
    init_sw_fields();
    sw_term((const su3**)hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
  }
  // we measure before the trajectory!
  if((mnl->rec_ev != 0) && (hf->traj_counter%mnl->rec_ev == 0)) {
    if(mnl->type != NDCLOVERRAT) phmc_compute_ev(hf->traj_counter-1, id, &Qtm_pm_ndbipsi);
    else phmc_compute_ev(hf->traj_counter-1, id, &Qsw_pm_ndbipsi);
  }

  // the Gaussian distributed random fields
  mnl->energy0 = 0.;
  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(mnl->pf, VOLUME/2, 1);

  random_spinor_field_eo(mnl->pf2, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 += square_norm(mnl->pf2, VOLUME/2, 1);

  solver_pm.max_iter = mnl->maxiter;
  solver_pm.squared_solver_prec = mnl->accprec;
  solver_pm.no_shifts = mnl->rat.np;
  solver_pm.shifts = mnl->rat.mu;
  solver_pm.type = CGMMSND;
  solver_pm.M_ndpsi = &Qtm_pm_ndpsi;
  if(mnl->type == NDCLOVERRATCOR) solver_pm.M_ndpsi = &Qsw_pm_ndpsi;
  solver_pm.sdim = VOLUME/2;
  solver_pm.rel_prec = g_relative_precision_flag;

  // apply B to the random field to generate pseudo-fermion fields
  assign(mnl->w_fields[0], mnl->pf, VOLUME/2);
  assign(mnl->w_fields[1], mnl->pf2, VOLUME/2);
  up0 = mnl->w_fields[0]; dn0 = mnl->w_fields[1];
  up1 = mnl->w_fields[2]; dn1 = mnl->w_fields[3];
	 
  for(int i = 1; i < 8; i++) {
    delta = apply_Z_ndpsi(up1, dn1, up0, dn0, id, hf, &solver_pm);
    assign_add_mul_r(mnl->pf, up1, coefs[i-1], VOLUME/2);
    assign_add_mul_r(mnl->pf2, dn1, coefs[i-1], VOLUME/2);
    if(delta < mnl->accprec) break;
    tup = up0; tdn = dn0;
    up0 = up1; dn0 = dn1;
    up1 = tup; dn1 = tdn;
  }
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) { 
      printf("called ndratcor_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}


double ndratcor_acc(const int id, hamiltonian_field_t * const hf) {
  solver_pm_t solver_pm;
  monomial * mnl = &monomial_list[id];
  double atime, etime, delta;
  spinor * up0, * dn0, * up1, * dn1, * tup, * tdn;
  double coefs[6] = {-1./2., 3./8., -5./16., 35./128., -63./256., 231./1024.};
  atime = gettime();
  nd_set_global_parameter(mnl);
  g_mu3 = 0.;
  if(mnl->type == NDCLOVERRATCOR) {
    sw_term((const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
  }
  mnl->energy1 = 0.;

  solver_pm.max_iter = mnl->maxiter;
  solver_pm.squared_solver_prec = mnl->accprec;
  solver_pm.no_shifts = mnl->rat.np;
  solver_pm.shifts = mnl->rat.mu;
  solver_pm.type = CGMMSND;
  solver_pm.M_ndpsi = &Qtm_pm_ndpsi;
  if(mnl->type == NDCLOVERRATCOR) solver_pm.M_ndpsi = &Qsw_pm_ndpsi;
  solver_pm.sdim = VOLUME/2;
  solver_pm.rel_prec = g_relative_precision_flag;

  // apply (Q R)^(-1) to pseudo-fermion fields
  assign(mnl->w_fields[4], mnl->pf, VOLUME/2);
  assign(mnl->w_fields[5], mnl->pf2, VOLUME/2);
  up0 = mnl->w_fields[0]; dn0 = mnl->w_fields[1];
  up1 = mnl->w_fields[2]; dn1 = mnl->w_fields[3];

  delta = apply_Z_ndpsi(up0, dn0, mnl->pf, mnl->pf2, id, hf, &solver_pm);
  assign_add_mul_r(mnl->w_fields[4], up0, coefs[0], VOLUME/2);
  assign_add_mul_r(mnl->w_fields[5], dn0, coefs[0], VOLUME/2);

  for(int i = 2; i < 8; i++) {
    if(delta < mnl->accprec) break;
    delta = apply_Z_ndpsi(up1, dn1, up0, dn0, id, hf, &solver_pm);
    assign_add_mul_r(mnl->w_fields[4], up1, coefs[i-1], VOLUME/2);
    assign_add_mul_r(mnl->w_fields[5], dn1, coefs[i-1], VOLUME/2);
    tup = up0; tdn = dn0;
    up0 = up1; dn0 = dn1;
    up1 = tup; dn1 = tdn;
  }

  mnl->energy1 = scalar_prod_r(mnl->pf, mnl->w_fields[4], VOLUME/2, 1);
  mnl->energy1 += scalar_prod_r(mnl->pf2, mnl->w_fields[5], VOLUME/2, 1);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) { // shoud be 3
      printf("called ndratcor_acc for id %d dH = %1.10e\n", id, mnl->energy1 - mnl->energy0);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}

// applies ((Q_h\tau_1 * R)^2 - 1)

double apply_Z_ndpsi(spinor * const k_up, spinor * const k_dn,
		     spinor * const l_up, spinor * const l_dn,
		     const int id, hamiltonian_field_t * const hf,
		     solver_pm_t * solver_pm) {
  monomial * mnl = &monomial_list[id];

  mnl->iter0 += cg_mms_tm_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
			     l_up, l_dn,
			     solver_pm);  
  
  // apply R to the pseudo-fermion fields
  assign(k_up, l_up, VOLUME/2);
  assign(k_dn, l_dn, VOLUME/2);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
    assign_add_mul_r(k_dn, g_chi_dn_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }

  // apply R a second time
  cg_mms_tm_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn,
	       solver_pm);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
    assign_add_mul_r(k_dn, g_chi_dn_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  mul_r(g_chi_up_spinor_field[mnl->rat.np], mnl->rat.A*mnl->rat.A, 
	k_up, VOLUME/2);
  mul_r(g_chi_dn_spinor_field[mnl->rat.np], mnl->rat.A*mnl->rat.A, 
	k_dn, VOLUME/2);
  // apply Q^2 and compute the residue
  solver_pm->M_ndpsi(k_up, k_dn,
		     g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np]);
  diff(k_up, k_up, l_up, VOLUME/2);
  diff(k_dn, k_dn, l_dn, VOLUME/2);
  double resi = square_norm(k_up, VOLUME/2, 1) + square_norm(k_dn, VOLUME/2, 1);
  if(g_debug_level > 2 && g_proc_id == 0) {
    printf("# NDRATCOR: ||Z * phi|| = %e\n", resi);
  }
  return(resi);
}

// computes ||(1 - C^dagger R C) phi||

void check_C_ndpsi(spinor * const k_up, spinor * const k_dn,
		   spinor * const l_up, spinor * const l_dn,
		   const int id, hamiltonian_field_t * const hf,
		   solver_pm_t * solver_pm) {
  monomial * mnl = &monomial_list[id];
  mnl->iter0 = cg_mms_tm_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
			     l_up, l_dn, solver_pm);

  assign(k_up, l_up, VOLUME/2);
  assign(k_dn, l_dn, VOLUME/2);

  // apply C to the random field to generate pseudo-fermion fields
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    // Q_h * tau^1 - i nu_j
    // this needs phmc_Cpol = 1 to work!
    if(mnl->type == NDCLOVERRATCOR || mnl->type == NDCLOVERRAT) {
      Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
			       g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			       I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    else {
      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			     I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    assign_add_mul(k_up, g_chi_up_spinor_field[mnl->rat.np], I*mnl->rat.rnu[j], VOLUME/2);
    assign_add_mul(k_dn, g_chi_dn_spinor_field[mnl->rat.np], I*mnl->rat.rnu[j], VOLUME/2);
  }
  //apply R
  solver_pm->shifts = mnl->rat.mu;
  cg_mms_tm_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn,
	       solver_pm);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
    assign_add_mul_r(k_dn, g_chi_dn_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  // apply C^dagger
  solver_pm->shifts = mnl->rat.nu;
  cg_mms_tm_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn, solver_pm);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    // Q_h * tau^1 + i nu_j
    if(mnl->type == NDCLOVERRATCOR || mnl->type == NDCLOVERRAT) {
      Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			     -I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    else {
      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			     -I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    assign_add_mul(k_up, g_chi_up_spinor_field[mnl->rat.np], -I*mnl->rat.rnu[j], VOLUME/2);
    assign_add_mul(k_dn, g_chi_dn_spinor_field[mnl->rat.np], -I*mnl->rat.rnu[j], VOLUME/2);
  }
  diff(k_up, k_up, l_up, VOLUME/2);  
  diff(k_dn, k_dn, l_dn, VOLUME/2);  
  double resi = square_norm(k_up, VOLUME/2, 1);
  resi += square_norm(k_dn, VOLUME/2, 1);
  if(g_proc_id == 0) printf("|| (1-C^dagger R C)*phi|| = %e\n", resi);

  return;
}
