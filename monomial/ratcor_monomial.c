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
#include "ratcor_monomial.h"

// computes ||(1 - C^dagger R C) phi||
void check_C_psi(spinor * const k_up, spinor * const l_up, 
		 const int id, hamiltonian_field_t * const hf,
		 solver_params_t * solver_params);

// applies (Q^2 R^2 -1) phi
double apply_Z_psi(spinor * const k_up, spinor * const l_up, 
		   const int id, hamiltonian_field_t * const hf,
		   solver_params_t * solver_params);



void ratcor_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime, delta;
  spinor * up0, * up1, * tup;
  double coefs[6] = {1./4., -3./32., 7./128., -77./2048., 231./8192., -1463./65536.}; // series of (1+x)^(1/4)
  double coefs_check[6] = {1./2., -1./8., 1./16., -5./128., 7./256., -21./1024.}; // series of (1+x)^(1/2)
  atime = gettime();
  nd_set_global_parameter(mnl);
  g_mu = 0.;
  g_mu3 = 0.;
  g_kappa = mnl->kappa;
  mnl->iter0 = 0;
  boundary(mnl->kappa);
  if(mnl->type == CLOVERRATCOR) {
    g_c_sw = mnl->c_sw;
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

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.type = mnl->solver;
  mnl->solver_params.M_psi = mnl->Qsq;
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;

  // apply B to the random field to generate pseudo-fermion fields
  assign(mnl->w_fields[0], mnl->pf, VOLUME/2);
  up0 = mnl->w_fields[0]; 
  up1 = mnl->w_fields[2]; 
	 
  for(int i = 1; i < 8; i++) {
    delta = apply_Z_psi(up1, up0, id, hf, &(mnl->solver_params) );
    assign_add_mul_r(mnl->pf, up1, coefs[i-1], VOLUME/2);
    if(delta < mnl->accprec) break;
    tup = up0;
    up0 = up1;
    up1 = tup;
  }
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called ratcor_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}


double ratcor_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime, delta;
  spinor * up0, * up1, * tup;
  double coefs[6] = {-1./2., 3./8., -5./16., 35./128., -63./256., 231./1024.};
  atime = gettime();
  nd_set_global_parameter(mnl);
  g_mu = 0.;
  g_mu3 = 0.;
  g_kappa = mnl->kappa;
  boundary(mnl->kappa);
  if(mnl->type == CLOVERRATCOR) {
    g_c_sw = mnl->c_sw;
    sw_term((const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert(EE, mnl->mu);
  }
  mnl->energy1 = 0.;

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.type = CGMMS;
  mnl->solver_params.M_psi = mnl->Qsq;
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;

  // apply (Q R)^(-1) to pseudo-fermion fields
  assign(mnl->w_fields[4], mnl->pf, VOLUME/2);
  up0 = mnl->w_fields[0];
  up1 = mnl->w_fields[2];

  delta = apply_Z_psi(up0, mnl->pf, id, hf, &(mnl->solver_params) );
  assign_add_mul_r(mnl->w_fields[4], up0, coefs[0], VOLUME/2);

  for(int i = 2; i < 8; i++) {
    if(delta < mnl->accprec) break;
    delta = apply_Z_psi(up1, up0, id, hf, &(mnl->solver_params) );
    assign_add_mul_r(mnl->w_fields[4], up1, coefs[i-1], VOLUME/2);
    tup = up0;
    up0 = up1;
    up1 = tup;
  }

  mnl->energy1 = scalar_prod_r(mnl->pf, mnl->w_fields[4], VOLUME/2, 1);
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called ratcor_acc for id %d dH = %1.10e\n", id, mnl->energy1 - mnl->energy0);
    }
  }
  return(mnl->energy1 - mnl->energy0);
}

// applies ((Q_h\tau_1 * R)^2 - 1)

double apply_Z_psi(spinor * const k_up,	spinor * const l_up, 
		     const int id, hamiltonian_field_t * const hf,
		     solver_params_t * solver_params) {
  monomial * mnl = &monomial_list[id];

  mnl->iter0 += solve_mms_tm(g_chi_up_spinor_field, l_up,
                             solver_params);  
  
  // apply R to the pseudo-fermion fields
  assign(k_up, l_up, VOLUME/2);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }

  // apply R a second time
  mnl->iter0 += solve_mms_tm(g_chi_up_spinor_field, k_up,
                             solver_params);

  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  mul_r(g_chi_up_spinor_field[mnl->rat.np], mnl->rat.A*mnl->rat.A, 
	k_up, VOLUME/2);

  // apply Q^2 
  solver_params->M_psi(k_up, g_chi_up_spinor_field[mnl->rat.np]);
  // compute the residue
  diff(k_up, k_up, l_up, VOLUME/2);
  double resi = square_norm(k_up, VOLUME/2, 1);
  
  if(g_debug_level > 2 && g_proc_id == 0) {
    printf("# RATCOR: ||Z * phi|| = %e\n", resi);
  }
  return(resi);
}

// computes ||(1 - C^dagger R C) phi||

void check_C_psi(spinor * const k_up, spinor * const l_up,
		 const int id, hamiltonian_field_t * const hf,
		 solver_params_t * solver_params) {
  monomial * mnl = &monomial_list[id];

  mnl->iter0 = solve_mms_tm(g_chi_up_spinor_field, l_up, solver_params);

  assign(k_up, l_up, VOLUME/2);

  // apply C to the random field to generate pseudo-fermion fields
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    if(mnl->type == CLOVERRATCOR || mnl->type == CLOVERRAT) {
      //Qsw_plus_psi(g_chi_up_spinor_field[mnl->rat.np], g_chi_up_spinor_field[j], 
      //			       I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    else {
      //Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
      //			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
      //		     I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    assign_add_mul(k_up, g_chi_up_spinor_field[mnl->rat.np], I*mnl->rat.rnu[j], VOLUME/2);
  }
  //apply R
  solver_params->shifts = mnl->rat.mu;
  mnl->iter0 += solve_mms_tm(g_chi_up_spinor_field, k_up,
                             solver_params);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  // apply C^dagger
  solver_params->shifts = mnl->rat.nu;
  mnl->iter0 += solve_mms_tm(g_chi_up_spinor_field, k_up,
	    solver_params);

  for(int j = (mnl->rat.np-1); j > -1; j--) {
    if(mnl->type == NDCLOVERRATCOR || mnl->type == NDCLOVERRAT) {
      //Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
      //			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
      //			     -I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    else {
      //Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np],
      //			     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
      //			     -I*mnl->rat.nu[j], 1., mnl->EVMaxInv);
    }
    assign_add_mul(k_up, g_chi_up_spinor_field[mnl->rat.np], -I*mnl->rat.rnu[j], VOLUME/2);
  }
  diff(k_up, k_up, l_up, VOLUME/2);  
  double resi = square_norm(k_up, VOLUME/2, 1);
  if(g_proc_id == 0) printf("|| (1-C^dagger R C)*phi|| = %e\n", resi);

  return;
}
