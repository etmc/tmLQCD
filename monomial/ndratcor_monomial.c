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
#include "solver/monomial_solve.h"
#include "deriv_Sb.h"
#include "init/init_chi_spinor_field.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/tm_operators_nd_32.h"
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
		   solver_params_t * solver_params);

// applies (Q^2 R^2 -1) phi
void apply_Z_ndpsi(spinor * const k_up, spinor * const k_dn,
		     spinor * const l_up, spinor * const l_dn,
		     const int id, hamiltonian_field_t * const hf,
		     solver_params_t * solver_params);



void ndratcor_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime, delta;
  spinor * up0, * dn0, * up1, * dn1, * tup, * tdn, * Zup, * Zdn;
  double coefs[6] = {1./4., -3./32., 7./128., -77./2048., 231./8192., -1463./65536.}; // series of (1+x)^(1/4)
  double coefs_check[6] = {1./2., -1./8., 1./16., -5./128., 7./256., -21./1024.}; // series of (1+x)^(1/2)
  atime = gettime();
  nd_set_global_parameter(mnl);
  g_mu3 = 0.;
  mnl->iter0 = 0;
  if(mnl->type == NDCLOVERRATCOR) {
    init_sw_fields();
    sw_term((const su3**)hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
    copy_32_sw_fields();
  }
  // we measure before the trajectory!
  if((mnl->rec_ev != 0) && (hf->traj_counter%mnl->rec_ev == 0)) {
    if(mnl->type != NDCLOVERRATCOR) phmc_compute_ev(hf->traj_counter-1, id, &Qtm_pm_ndbipsi);
    else phmc_compute_ev(hf->traj_counter-1, id, &Qsw_pm_ndbipsi);
  }

  // the Gaussian distributed random fields
  mnl->energy0 = 0.;
  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(mnl->pf, VOLUME/2, 1);

  random_spinor_field_eo(mnl->pf2, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 += square_norm(mnl->pf2, VOLUME/2, 1);

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.type = mnl->solver;
  mnl->solver_params.M_ndpsi = &Qtm_pm_ndpsi;
  mnl->solver_params.M_ndpsi32 = &Qtm_pm_ndpsi_32;    
  if(mnl->type == NDCLOVERRATCOR) {
    mnl->solver_params.M_ndpsi = &Qsw_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = &Qsw_pm_ndpsi_32;
  }
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;

  // apply B to the random field to generate pseudo-fermion fields
  up0 = mnl->w_fields[0]; dn0 = mnl->w_fields[1];
  up1 = mnl->w_fields[2]; dn1 = mnl->w_fields[3];
  Zup = mnl->w_fields[4]; Zdn = mnl->w_fields[5];

  apply_Z_ndpsi(up0, dn0, mnl->pf, mnl->pf2, id, hf, &(mnl->solver_params));
  // computing correction to energy1
  delta = coefs_check[0]*(scalar_prod_r(mnl->pf, up0, VOLUME/2, 1) + scalar_prod_r(mnl->pf2, dn0, VOLUME/2, 1));
  if(g_debug_level > 2 && g_proc_id == 0)
    printf("# NDRATCOR heatbath: c_%d*(R * Z^%d * R) = %e\n", 1, 1, delta);
  // debug for showing that the old check was giving a smaller delta
  if(g_debug_level > 3) {
    double delta_old = square_norm(up0, VOLUME/2, 1) + square_norm(dn0, VOLUME/2, 1);
    if(g_proc_id == 0) {
      printf("# NDRATCOR old check: || Z^%d * R ||^2 = %e\n", 1, delta_old);
      printf("# NDRATCOR new check: (c_%d*(R * Z^%d * R))^2 = %e\n", 1, 1, delta*delta);
    }
  }

  if(delta*delta > mnl->accprec) {
    assign_add_mul_r(mnl->pf, up0, coefs[0], VOLUME/2);
    assign_add_mul_r(mnl->pf2, dn0, coefs[0], VOLUME/2);
    
    // saving first application
    assign(Zup, up0, VOLUME/2);
    assign(Zdn, dn0, VOLUME/2);
    
    
    for(int i = 2; i < 8; i++) {
      // computing next order correction to energy1
      delta = coefs_check[i-1]*(scalar_prod_r(Zup, up0, VOLUME/2, 1) + scalar_prod_r(Zup, dn0, VOLUME/2, 1)); 
      if(g_debug_level > 2 && g_proc_id == 0)
        printf("# NDRATCOR heatbath: c_%d*(R * Z^%d * R) = %e\n", i, i, delta);
      // debug for showing that the old check was giving a smaller delta
      if(g_debug_level > 3) {
        double delta_old = square_norm(up0, VOLUME/2, 1) + square_norm(dn0, VOLUME/2, 1);
        if(g_proc_id == 0) {
          printf("# NDRATCOR old check: || Z^%d * R ||^2 = %e\n", 1, delta_old);
          printf("# NDRATCOR new check: (c_%d*(R * Z^%d * R))^2 = %e\n", 1, 1, delta*delta);
        }
      }
      if(delta*delta < mnl->accprec) break;

      apply_Z_ndpsi(up1, dn1, up0, dn0, id, hf, &(mnl->solver_params));
      
      assign_add_mul_r(mnl->pf, up1, coefs[i-1], VOLUME/2);
      assign_add_mul_r(mnl->pf2, dn1, coefs[i-1], VOLUME/2);

      tup = up0; tdn = dn0;
      up0 = up1; dn0 = dn1;
      up1 = tup; dn1 = tdn;
    }
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
    copy_32_sw_fields();
  }
  mnl->energy1 = square_norm(mnl->pf, VOLUME/2, 1) + square_norm(mnl->pf2, VOLUME/2, 1);

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.type = mnl->solver;
  mnl->solver_params.M_ndpsi = &Qtm_pm_ndpsi;
  mnl->solver_params.M_ndpsi32 = &Qtm_pm_ndpsi_32;    
  if(mnl->type == NDCLOVERRATCOR) {
    mnl->solver_params.M_ndpsi = &Qsw_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = &Qsw_pm_ndpsi_32;
  }
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;

  // apply (Q R)^(-1) to pseudo-fermion fields
  up0 = mnl->w_fields[0]; dn0 = mnl->w_fields[1];
  up1 = mnl->w_fields[2]; dn1 = mnl->w_fields[3];

  apply_Z_ndpsi(up0, dn0, mnl->pf, mnl->pf2, id, hf, &(mnl->solver_params));
  delta = coefs[0]*(scalar_prod_r(mnl->pf, up0, VOLUME/2, 1) + scalar_prod_r(mnl->pf2, dn0, VOLUME/2, 1));
  mnl->energy1 += delta;
  if(g_debug_level > 2 && g_proc_id == 0)
    printf("# NDRATCOR acc step: c_%d*(phi * Z^%d * phi) = %e\n", 1, 1, delta);

  for(int i = 2; i < 8; i++) {
    if(delta*delta < mnl->accprec) break;

    delta = coefs[i-1]*(square_norm(up0, VOLUME/2, 1) + square_norm(dn0, VOLUME/2, 1)); 
    mnl->energy1 += delta;
    if(g_debug_level > 2 && g_proc_id == 0)
      printf("# NDRATCOR acc step: c_%d*(phi * Z^%d * phi) = %e\n", i, i, delta);
    i++; //incrementing i
    if(delta*delta < mnl->accprec) break;

    apply_Z_ndpsi(up1, dn1, up0, dn0, id, hf, &(mnl->solver_params));
    delta = coefs[i-1]*(scalar_prod_r(up0, up1, VOLUME/2, 1) + scalar_prod_r(dn0, dn1, VOLUME/2, 1));
    mnl->energy1 += delta;
    if(g_debug_level > 2 && g_proc_id == 0)
      printf("# NDRATCOR acc step: c_%d*(phi * Z^%d * phi) = %e\n", i, i, delta);

    tup = up0; tdn = dn0;
    up0 = up1; dn0 = dn1;
    up1 = tup; dn1 = tdn;
  }

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
void apply_Z_ndpsi(spinor * const k_up, spinor * const k_dn,
		     spinor * const l_up, spinor * const l_dn,
		     const int id, hamiltonian_field_t * const hf,
		     solver_params_t * solver_params) {
  monomial * mnl = &monomial_list[id];

  mnl->iter0 += solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
			                       l_up, l_dn, solver_params);  
  
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
  mnl->iter0 += solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn,
	       solver_params);
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
  solver_params->M_ndpsi(k_up, k_dn,
		     g_chi_up_spinor_field[mnl->rat.np], g_chi_dn_spinor_field[mnl->rat.np]);
  diff(k_up, k_up, l_up, VOLUME/2);
  diff(k_dn, k_dn, l_dn, VOLUME/2);
  
}

// computes ||(1 - C^dagger R C) phi||
void check_C_ndpsi(spinor * const k_up, spinor * const k_dn,
		   spinor * const l_up, spinor * const l_dn,
		   const int id, hamiltonian_field_t * const hf,
		   solver_params_t * solver_params) {
  monomial * mnl = &monomial_list[id];
  mnl->iter0 = solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
			     l_up, l_dn, solver_params);

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
  solver_params->shifts = mnl->rat.mu;
  solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn,
	       solver_params);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(k_up, g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
    assign_add_mul_r(k_dn, g_chi_dn_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  // apply C^dagger
  solver_params->shifts = mnl->rat.nu;
  solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
	       k_up, k_dn, solver_params);
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
