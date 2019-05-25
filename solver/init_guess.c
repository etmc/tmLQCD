/***********************************************************************
 *
 *
 * Copyright (C) 2016 Simone Bacchio
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
#include "global.h"
#include "su3.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/solver.h"
#include "solver_field.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_nd.h"
#include "init_guess.h"
#include <io/params.h>

int init_guess_mms(spinor ** const P, spinor * const Q,
                   int shift, solver_params_t * const solver_params) {
  double * shifts=solver_params->shifts;
  int no_shifts = solver_params->no_shifts;
  if(shift==no_shifts-1) {
    zero_spinor_field(P[shift], solver_params->sdim);
  } else {
    double coeff;
    for( int j = no_shifts-1; j > shift; j-- ) {
      coeff = 1;
      for( int k = no_shifts-1; k > shift; k-- ) {
        if(j!=k)
          coeff *= (shifts[k]*shifts[k]-shifts[shift]*shifts[shift])/
                   (shifts[k]*shifts[k]-shifts[j]*shifts[j]);
      }
      if(j==no_shifts-1) {
        mul_r(P[shift], coeff, P[j], solver_params->sdim);
      } else {
        assign_add_mul_r(P[shift], P[j], coeff, solver_params->sdim);
      }
    }
  }
  if(g_debug_level > 2){
    double old_g_mu3 = g_mu3;
    spinor** temp;
    if(solver_params->sdim == VOLUME/2) {
      init_solver_field(&temp, VOLUMEPLUSRAND/2, 1);
    } else {
      init_solver_field(&temp, VOLUMEPLUSRAND, 1);
    }

    g_mu3 = solver_params->shifts[shift]; 
    solver_params->M_psi( temp[0], P[shift]);
    g_mu3 = old_g_mu3;

    diff( temp[0], temp[0], Q, solver_params->sdim);
    double res = sqrt(square_norm(temp[0], solver_params->sdim, 1)/square_norm(Q, solver_params->sdim, 1));
      
    finalize_solver(temp, 1);
    if(g_proc_id == 0)
      printf("INITIAL GUESS: shift id=%d value=%e  relative residual: %e\n",shift,shifts[shift],res); 
  }

}

int init_guess_mms_nd(spinor ** const Pup, spinor ** const Pdn, 
                      spinor * const Qup, spinor * const Qdn, 
                      int shift, solver_params_t * solver_params) {
  double * shifts=solver_params->shifts;
  int no_shifts = solver_params->no_shifts;
  if(shift==no_shifts-1) {
    zero_spinor_field(Pup[shift], solver_params->sdim);
    zero_spinor_field(Pdn[shift], solver_params->sdim);
  } else {
    double coeff;
    for( int j = no_shifts-1; j > shift; j-- ) {
      coeff = 1;
      for( int k = no_shifts-1; k > shift; k-- ) {
        if(j!=k)
          coeff *= (shifts[k]*shifts[k]-shifts[shift]*shifts[shift])/
                   (shifts[k]*shifts[k]-shifts[j]*shifts[j]);
      }
      if(j==no_shifts-1) {
        mul_r(Pup[shift], coeff, Pup[j], solver_params->sdim);
        mul_r(Pdn[shift], coeff, Pdn[j], solver_params->sdim);
      } else {
        assign_add_mul_r(Pup[shift], Pup[j], coeff, solver_params->sdim);
        assign_add_mul_r(Pdn[shift], Pdn[j], coeff, solver_params->sdim);
      }
    }
  }
  if(g_debug_level > 2){
    double old_g_shift = g_shift;
    matrix_mult_nd f = Qtm_pm_ndpsi_shift;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi ) 
      f = Qsw_pm_ndpsi_shift;
    spinor** temp;
    if(solver_params->sdim == VOLUME/2) {
      init_solver_field(&temp, VOLUMEPLUSRAND/2, 2);
    } else {
      init_solver_field(&temp, VOLUMEPLUSRAND, 2);
    }

    g_shift = shifts[shift]*shifts[shift]; 
    f( temp[0], temp[1], Pup[shift], Pdn[shift]);
    g_shift = old_g_shift;

    diff( temp[0], temp[0], Qup, solver_params->sdim);
    diff( temp[1], temp[1], Qdn, solver_params->sdim);
    double res = sqrt(square_norm(temp[0], solver_params->sdim, 1)+square_norm(temp[1], solver_params->sdim, 1))/
      sqrt(square_norm(Qup, solver_params->sdim, 1)+square_norm(Qdn, solver_params->sdim, 1));
      
    finalize_solver(temp, 2);
    if(g_proc_id == 0)
      printf("INITIAL GUESS ND: shift id=%d value=%e  relative residual: %e\n",shift,shifts[shift],res); 
  }
}

int init_guess_mms_nd_plus(spinor ** const Pup, spinor ** const Pdn, 
                           spinor * const Qup, spinor * const Qdn, 
                           int shift, solver_params_t * solver_params) {
  double * shifts=solver_params->shifts;
  int no_shifts = solver_params->no_shifts;
  if(shift==no_shifts-1) {
    zero_spinor_field(Pup[shift], solver_params->sdim);
    zero_spinor_field(Pdn[shift], solver_params->sdim);
  } else {
    double coeff;
    for( int j = no_shifts-1; j > shift; j-- ) {
      coeff = 1;
      for( int k = no_shifts-1; k > shift; k-- ) {
        if(j!=k)
          coeff *= (shifts[k]-shifts[shift])/(shifts[k]-shifts[j]);
      }
      if(j==no_shifts-1) {
        mul_r(Pup[shift], coeff, Pup[j], solver_params->sdim);
        mul_r(Pdn[shift], coeff, Pdn[j], solver_params->sdim);
      } else {
        assign_add_mul_r(Pup[shift], Pup[j], coeff, solver_params->sdim);
        assign_add_mul_r(Pdn[shift], Pdn[j], coeff, solver_params->sdim);
      }
    }
  }
  if(g_debug_level > 2){
    double old_g_shift = g_shift;
    matrix_mult_nd f = Qtm_tau1_ndpsi_add_Ishift;
    if( solver_params->M_ndpsi == Qsw_pm_ndpsi )
      f = Qsw_tau1_ndpsi_add_Ishift;
    spinor** temp;
    if(solver_params->sdim == VOLUME/2) {
      init_solver_field(&temp, VOLUMEPLUSRAND/2, 2);
    } else {
      init_solver_field(&temp, VOLUMEPLUSRAND, 2);
    }

    g_shift = shifts[shift]*shifts[shift]; 
    f( temp[0], temp[1], Pup[shift], Pdn[shift]);
    g_shift = old_g_shift;

    diff( temp[0], temp[0], Qup, solver_params->sdim);
    diff( temp[1], temp[1], Qdn, solver_params->sdim);
    double res = sqrt(square_norm(temp[0], solver_params->sdim, 1)+square_norm(temp[1], solver_params->sdim, 1))/
      sqrt(square_norm(Qup, solver_params->sdim, 1)+square_norm(Qdn, solver_params->sdim, 1));
      
    finalize_solver(temp, 2);
    if(g_proc_id == 0)
      printf("INITIAL GUESS ND PLUS: shift id=%d value=%e  relative residual: %e\n",shift,shifts[shift],res); 
  }
}
