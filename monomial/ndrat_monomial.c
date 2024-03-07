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
# include<tmlqcd_config.h>
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
#include "default_input_values.h"

void nd_set_global_parameter(monomial * const mnl) {

  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  g_kappa = mnl->kappa;
  g_c_sw = mnl->c_sw;
  boundary(g_kappa);
  phmc_cheb_evmin = mnl->EVMin;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmax = mnl->EVMax;
  phmc_Cpol = 1.;
  // used for preconditioning in cloverdetrat
  g_mu3 = 0.;

  return;
}



void print_spinor_updow_shift( double** up, double** down,   int odd, int num_shifts) {
  tm_stopwatch_push(&g_timers, __func__, "");

  const int change_sign[4] = {-1, 1, 1, -1};
  const int change_spin[4] = {3, 2, 1, 0};
  const int Vh = VOLUME/2;

  // now copy and reorder from tempSpinor to spinor
  for(int shift = 0; shift < num_shifts; shift++){
  for(int ud = 0; ud < 2; ud++){
    double **sp;
    if (ud==0) sp=up;
    if (ud==1) sp=down;
    for( int x0=0; x0<T; x0++ )
      for( int x1=0; x1<LX; x1++ )
        for( int x2=0; x2<LY; x2++ )
          for( int x3=0; x3<LZ; x3++ ) {
            const int q_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
            const int tm_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;

            const int oddBit = (x0+x1+x2+x3) & 1;
            if( oddBit == odd ){
              for(int tm_spin = 0; tm_spin < 4; tm_spin++){
                for(int col = 0; col < 3; col++){
                  if (g_proc_id==0) printf("MARCOfrom TMLQCD (%d %d %d %d),  %d %d, s=%d  ud=%d    %g  %g\n",x0,x1,x2,x3, tm_spin, col ,
                    shift, ud,
                    sp[shift][24*tm_eo_idx + 6*tm_spin + 2*col + 0],
                    sp[shift][24*tm_eo_idx + 6*tm_spin + 2*col + 1]);
                }
              }
            }
          }
  }
  }
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}      

void print_spinor_updow( double* up, double* down,   int odd , int shift) {
  tm_stopwatch_push(&g_timers, __func__, "");

  const int change_sign[4] = {-1, 1, 1, -1};
  const int change_spin[4] = {3, 2, 1, 0};
  const int Vh = VOLUME/2;

  // now copy and reorder from tempSpinor to spinor
  
  for(int ud = 0; ud < 2; ud++){
    double *sp;
    if (ud==0) sp=up;
    if (ud==1) sp=down;
    for( int x0=0; x0<T; x0++ )
      for( int x1=0; x1<LX; x1++ )
        for( int x2=0; x2<LY; x2++ )
          for( int x3=0; x3<LZ; x3++ ) {
            const int q_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
            const int tm_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;

            const int oddBit = (x0+x1+x2+x3) & 1;
            if( oddBit == odd ){
              for(int tm_spin = 0; tm_spin < 4; tm_spin++){
                for(int col = 0; col < 3; col++){
                  if (g_proc_id==0) printf("MARCOfrom TMLQCD (%d %d %d %d),  %d %d, s=%d  ud=%d    %g  %g\n",x0,x1,x2,x3, tm_spin, col ,
                    shift, ud,
                    sp[24*tm_eo_idx + 6*tm_spin + 2*col + 0],
                    sp[24*tm_eo_idx + 6*tm_spin + 2*col + 1]);
                }
              }
            }
  
  }
  }
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}      


void set_spinor_updow_shift( double** up, double** down,   int odd, int num_shifts) {
  tm_stopwatch_push(&g_timers, __func__, "");

  const int change_sign[4] = {-1, 1, 1, -1};
  const int change_spin[4] = {3, 2, 1, 0};
  const int Vh = VOLUME/2;

  // now copy and reorder from tempSpinor to spinor
  for(int shift = 0; shift < num_shifts; shift++){
  for(int ud = 0; ud < 2; ud++){
    double **sp;
    if (ud==0) sp=up;
    if (ud==1) sp=down;
    for( int x0=0; x0<T; x0++ )
      for( int x1=0; x1<LX; x1++ )
        for( int x2=0; x2<LY; x2++ )
          for( int x3=0; x3<LZ; x3++ ) {
            const int q_eo_idx = (x1 + LX*x2 + LY*LX*x3 + LZ*LY*LX*x0)/2;
            const int tm_eo_idx = (x3 + LZ*x2 + LY*LZ*x1 + LX*LY*LZ*x0)/2;

            const int oddBit = (x0+x1+x2+x3) & 1;
            if( oddBit == odd ){
              for(int q_spin = 0; q_spin < 4; q_spin++){
                const int tm_spin = change_spin[q_spin];
                for(int col = 0; col < 3; col++){
                    sp[shift][24*tm_eo_idx + 6*tm_spin + 2*col + 0]=0;
                    sp[shift][24*tm_eo_idx + 6*tm_spin + 2*col + 1 + Vh]=0;
                }
              }
            }
          }
  }
  }
  
  up[0][24*0 + 6*0 + 2*0 + 0]=1;
  up[0][24*0 + 6*0 + 2*0 + 1]=3;
  up[0][24*0 + 6*0 + 2*1 + 0]=4;
  down[0][24*0 + 6*0 + 2*0 + 0]=2;
  tm_stopwatch_pop(&g_timers, 0, 0, "TM_QUDA");
}      

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void ndrat_derivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  nd_set_global_parameter(mnl);
  if(mnl->type == NDCLOVERRAT) {
    tm_stopwatch_push(&g_timers, "su3_zero", "");
#ifdef TM_USE_OMP
    #pragma omp parallel for
#endif
    for(int i = 0; i < VOLUME; i++) { 
      for(int mu = 0; mu < 4; mu++) { 
        _su3_zero(swm[i][mu]);
        _su3_zero(swp[i][mu]);
      }
    }
    tm_stopwatch_pop(&g_timers, 0, 1, "");
  
    // we compute the clover term (1 + T_ee(oo)) for all sites x
    sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    // we invert it for the even sites only
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
    copy_32_sw_fields();
  }
  mnl->forcefactor = mnl->EVMaxInv;

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->forceprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.rel_prec = g_relative_precision_flag;
  mnl->solver_params.type = mnl->solver; 
  mnl->solver_params.M_ndpsi = &Qtm_pm_ndpsi;
  mnl->solver_params.M_ndpsi32 = &Qtm_pm_ndpsi_32;    
  if(mnl->type == NDCLOVERRAT) {
    mnl->solver_params.M_ndpsi = &Qsw_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = &Qsw_pm_ndpsi_32;
  }
  mnl->solver_params.sdim = VOLUME/2;

  // this generates all X_j,o (odd sites only) -> g_chi_up|dn_spinor_field
  mnl->iter1 += solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
                             mnl->pf, mnl->pf2, &(mnl->solver_params) );

  // set_spinor_updow_shift(  g_chi_up_spinor_field,  g_chi_dn_spinor_field,   1, mnl->rat.np);

  if ( mnl->external_library == QUDA_LIB){
    if(!mnl->even_odd_flag) {
      fatal_error("QUDA support only even_odd_flag",__func__);
    }
#ifdef TM_USE_QUDA
    if (g_debug_level > 3) {
#ifdef TM_USE_MPI
      // FIXME: here we do not need to set to zero the interior but only the halo
#ifdef TM_USE_OMP
      # pragma omp parallel for
#endif // end setting to zero the halo when using MPI
      for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand);i++) { 
        for(int mu=0;mu<4;mu++) { 
          _zero_su3adj(debug_derivative[i][mu]);
        }
      }
#endif
      // we copy only the interior
      memcpy(debug_derivative[0], hf->derivative[0], 4*VOLUME*sizeof(su3adj));
    }
    if(! (mnl->type == NDCLOVERRAT)){
      fatal_error("QUDA support only mnl->type = NDCLOVERRAT",__func__);
    }
    printf("are the same? mnl->rat.np =%d  mnl->solver_params->no_shifts=%d\n",mnl->rat.np,  mnl->solver_params.no_shifts);
    compute_ndcloverrat_derivative_quda(mnl, hf, g_chi_up_spinor_field, g_chi_dn_spinor_field, &(mnl->solver_params) );

    if (g_debug_level > 3){
      su3adj **given = hf->derivative;
      hf->derivative = debug_derivative;
      // int store_solver = mnl->solver;
      // mnl->solver_params.external_inverter = NO_EXT_INV;
      // mnl->solver = CGMMSND;
      mnl->external_library = NO_EXT_LIB;
      tm_debug_printf( 0, 3, "Recomputing the derivative from tmLQCD\n");
      ndrat_derivative(id, hf);
      #ifdef TM_USE_MPI
        xchange_deri(hf->derivative);// this function use ddummy inside
      #endif
      compare_derivative(mnl, given, hf->derivative, 1e-9, "ndrat_derivative");
      mnl->external_library = QUDA_LIB;
      hf->derivative = given;
    }
  #endif // no other option, TM_USE_QUDA already checked by solver
  }
  else{
    print_spinor_updow_shift( g_chi_up_spinor_field, g_chi_dn_spinor_field,   1, mnl->rat.np);
    for(int j = 0; j < mnl->rat.np; j++) {
      if(mnl->type == NDCLOVERRAT) {
        // multiply with Q_h * tau^1 + i mu_j to get Y_j,o (odd sites)
        // needs phmc_Cpol = 1 to work for ndrat!
        tm_stopwatch_push(&g_timers, "Qsw_tau1_sub_const_ndpsi", "");
        Qsw_tau1_sub_const_ndpsi(mnl->w_fields[0], mnl->w_fields[1],
              g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
              -I*mnl->rat.mu[j], 1., mnl->EVMaxInv);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
        
        /* Get the even parts X_j,e */
        /* H_eo_... includes tau_1 */
        tm_stopwatch_push(&g_timers, "H_eo_sw_ndpsi", "");
        H_eo_sw_ndpsi(mnl->w_fields[2], mnl->w_fields[3], 
          g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j]);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
        // print_spinor_updow( mnl->w_fields[2], mnl->w_fields[3],   0, j);

      } else {
        // multiply with Q_h * tau^1 + i mu_j to get Y_j,o (odd sites)
        // needs phmc_Cpol = 1 to work for ndrat!
        tm_stopwatch_push(&g_timers, "Q_tau1_sub_const_ndpsi", "");
        Q_tau1_sub_const_ndpsi(mnl->w_fields[0], mnl->w_fields[1],
            g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
            -I*mnl->rat.mu[j], 1., mnl->EVMaxInv);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
        
        /* Get the even parts X_j,e */
        /* H_eo_... includes tau_1 */
        tm_stopwatch_push(&g_timers, "H_eo_tm_ndpsi", "");
        H_eo_tm_ndpsi(mnl->w_fields[2], mnl->w_fields[3], 
          g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], EO);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
      }
      /* X_j,e^dagger \delta M_eo Y_j,o */
      deriv_Sb(EO, mnl->w_fields[2], mnl->w_fields[0], 
        hf, mnl->rat.rmu[j]*mnl->forcefactor);
      deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1],
        hf, mnl->rat.rmu[j]*mnl->forcefactor);

      if(mnl->type == NDCLOVERRAT) {
        /* Get the even parts Y_j,e */
        tm_stopwatch_push(&g_timers, "H_eo_sw_ndpsi", "");
        H_eo_sw_ndpsi(mnl->w_fields[4], mnl->w_fields[5], mnl->w_fields[0], mnl->w_fields[1]);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
      }
      else {
        /* Get the even parts Y_j,e */
        tm_stopwatch_push(&g_timers, "H_eo_tm_ndpsi", "");
        H_eo_tm_ndpsi(mnl->w_fields[4], mnl->w_fields[5], mnl->w_fields[0], mnl->w_fields[1], EO);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
      }

      /* X_j,o \delta M_oe Y_j,e */
      deriv_Sb(OE, g_chi_up_spinor_field[j], mnl->w_fields[4], 
        hf, mnl->rat.rmu[j]*mnl->forcefactor);
      deriv_Sb(OE, g_chi_dn_spinor_field[j], mnl->w_fields[5], 
        hf, mnl->rat.rmu[j]*mnl->forcefactor);

      if(mnl->type == NDCLOVERRAT) {
        // even/even sites sandwiched by tau_1 gamma_5 Y_e and gamma_5 X_e
        sw_spinor_eo(EE, mnl->w_fields[5], mnl->w_fields[2], 
        mnl->rat.rmu[j]*mnl->forcefactor);
        // odd/odd sites sandwiched by tau_1 gamma_5 Y_o and gamma_5 X_o
        sw_spinor_eo(OO, g_chi_up_spinor_field[j], mnl->w_fields[1],
        mnl->rat.rmu[j]*mnl->forcefactor);
        
        // even/even sites sandwiched by tau_1 gamma_5 Y_e and gamma_5 X_e
        sw_spinor_eo(EE, mnl->w_fields[4], mnl->w_fields[3], 
        mnl->rat.rmu[j]*mnl->forcefactor);
        // odd/odd sites sandwiched by tau_1 gamma_5 Y_o and gamma_5 X_o
        sw_spinor_eo(OO, g_chi_dn_spinor_field[j], mnl->w_fields[0],
        mnl->rat.rmu[j]*mnl->forcefactor);
      }
    }
    // trlog part does not depend on the normalisation
    if(mnl->type == NDCLOVERRAT && mnl->trlog) {
      sw_deriv_nd(EE);
    }
    if(mnl->type == NDCLOVERRAT) {
      sw_all(hf, mnl->kappa, mnl->c_sw);
    }
  }
  
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}


void ndrat_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  nd_set_global_parameter(mnl);
  mnl->iter1 = 0;
  if(mnl->type == NDCLOVERRAT) {
    init_sw_fields();
    sw_term((const su3**)hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
    copy_32_sw_fields();
  }
  // we measure before the trajectory!
  if((mnl->rec_ev != 0) && (hf->traj_counter%mnl->rec_ev == 0)) {
    if(mnl->type != NDCLOVERRAT) phmc_compute_ev(hf->traj_counter-1, id, &Qtm_pm_ndbipsi);
    else phmc_compute_ev(hf->traj_counter-1, id, &Qsw_pm_ndbipsi);
  }

  // the Gaussian distributed random fields
  tm_stopwatch_push(&g_timers, "random_energy0", "");
  mnl->energy0 = 0.;
  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(mnl->pf, VOLUME/2, 1);

  random_spinor_field_eo(mnl->pf2, mnl->rngrepro, RN_GAUSS);
  mnl->energy0 += square_norm(mnl->pf2, VOLUME/2, 1);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  
  // set solver parameters
  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.nu;
  mnl->solver_params.type = mnl->solver;
  mnl->solver_params.M_ndpsi = &Qtm_pm_ndpsi;
  mnl->solver_params.M_ndpsi32 = &Qtm_pm_ndpsi_32;    
  if(mnl->type == NDCLOVERRAT) {
    mnl->solver_params.M_ndpsi = &Qsw_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = &Qsw_pm_ndpsi_32;
  }
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;
  mnl->iter0 = solve_mms_nd_plus(g_chi_up_spinor_field, g_chi_dn_spinor_field,
                                 mnl->pf, mnl->pf2, &(mnl->solver_params) );

  assign(mnl->w_fields[2], mnl->pf, VOLUME/2);
  assign(mnl->w_fields[3], mnl->pf2, VOLUME/2);

  // apply C to the random field to generate pseudo-fermion fields
  for(int j = (mnl->rat.np-1); j > -1; j--) {
      assign_add_mul(mnl->pf, g_chi_up_spinor_field[j], I*mnl->rat.rnu[j], VOLUME/2);
      assign_add_mul(mnl->pf2, g_chi_dn_spinor_field[j], I*mnl->rat.rnu[j], VOLUME/2);
  }

  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called ndrat_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}


double ndrat_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  nd_set_global_parameter(mnl);
  if(mnl->type == NDCLOVERRAT) {
    sw_term((const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
    copy_32_sw_fields();
  }
  mnl->energy1 = 0.;

  mnl->solver_params.max_iter = mnl->maxiter;
  mnl->solver_params.squared_solver_prec = mnl->accprec;
  mnl->solver_params.no_shifts = mnl->rat.np;
  mnl->solver_params.shifts = mnl->rat.mu;
  mnl->solver_params.type = mnl->solver;
  
  mnl->solver_params.M_ndpsi = &Qtm_pm_ndpsi;
  mnl->solver_params.M_ndpsi32 = &Qtm_pm_ndpsi_32; 
  if(mnl->type == NDCLOVERRAT) {
    mnl->solver_params.M_ndpsi = &Qsw_pm_ndpsi;
    mnl->solver_params.M_ndpsi32 = &Qsw_pm_ndpsi_32;
  }
  mnl->solver_params.sdim = VOLUME/2;
  mnl->solver_params.rel_prec = g_relative_precision_flag;
  mnl->iter0 += solve_mms_nd(g_chi_up_spinor_field, g_chi_dn_spinor_field,
                            mnl->pf, mnl->pf2, &(mnl->solver_params) );

  // apply R to the pseudo-fermion fields
  assign(mnl->w_fields[0], mnl->pf, VOLUME/2);
  assign(mnl->w_fields[1], mnl->pf2, VOLUME/2);
  for(int j = (mnl->rat.np-1); j > -1; j--) {
    assign_add_mul_r(mnl->w_fields[0], g_chi_up_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
    assign_add_mul_r(mnl->w_fields[1], g_chi_dn_spinor_field[j], 
		     mnl->rat.rmu[j], VOLUME/2);
  }
  
  tm_stopwatch_push(&g_timers, "scalar_prod_r", "");
  mnl->energy1 = scalar_prod_r(mnl->pf, mnl->w_fields[0], VOLUME/2, 1);
  mnl->energy1 += scalar_prod_r(mnl->pf2, mnl->w_fields[1], VOLUME/2, 1);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called ndrat_acc for id %d, H_1 = %.10e, dH = %1.10e\n", id, mnl->energy1,  mnl->energy1 - mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return(mnl->energy1 - mnl->energy0);
}


int init_ndrat_monomial(const int id) {
  monomial * mnl = &monomial_list[id];  
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  int scale = 0;

  if(mnl->type == RAT || mnl->type == CLOVERRAT ||
     mnl->type == RATCOR || mnl->type == CLOVERRATCOR) 
    scale = 1;

  if(scale) {
    // When scale = 1 
    //   the rational approximation is done for the standard operator 
    //   which have eigenvalues between EVMin and EVMax.  Indeed the 
    //   parameters of the rational approximation are scaled. Thus 
    //   additional scaling of the operator (EVMaxInv) is not required.
    mnl->EVMin = mnl->StildeMin;
    mnl->EVMax = mnl->StildeMax;
    mnl->EVMaxInv = 1.;
  } else {
    // When scale = 0 
    //   the rational approximation is done for the normalized operator 
    //   which have eigenvalues between EVMin/EVMax and 1. Thus the 
    //   operator need to be scaled by EVMaxInv=1/EVMax.
    mnl->EVMin = mnl->StildeMin / mnl->StildeMax;
    mnl->EVMax = 1.;
    mnl->EVMaxInv = 1./sqrt(mnl->StildeMax);
  }

  init_rational(&mnl->rat, scale);

  if(mnl->type == RAT || mnl->type == CLOVERRAT ||
     mnl->type == RATCOR || mnl->type == CLOVERRATCOR) {
    if(init_chi_spinor_field(VOLUMEPLUSRAND/2, (mnl->rat.np+2)/2) != 0) {
      fprintf(stderr, "Not enough memory for Chi fields! Aborting...\n");
      exit(0);
    }
  } else {
    if(init_chi_spinor_field(VOLUMEPLUSRAND/2, (mnl->rat.np+1)) != 0) {
      fprintf(stderr, "Not enough memory for Chi fields! Aborting...\n");
      exit(0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return(0);
}

