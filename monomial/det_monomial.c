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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "deriv_Sb_D_psi.h"
#include "operator/tm_operators.h"
#include "operator/Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "solver/solver.h"
#include "solver/monomial_solve.h"
#include "read_input.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "monomial/monomial.h"
#include "det_monomial.h"

/* think about chronological solver ! */

void det_derivative(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  mnl->forcefactor = 1.;
  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);
  
  g_mu = mnl->mu;
  g_kappa = mnl->kappa;
  boundary(g_kappa);

  if(mnl->even_odd_flag) {
    /*********************************************************************
     * 
     * even/odd version 
     *
     * This a term is det(\hat Q^2(\mu))
     *
     *********************************************************************/
    
    /* Invert Q_{+} Q_{-} */
    /* X_o -> w_fields[1] */
    chrono_guess(mnl->w_fields[1], mnl->pf, mnl->csg_field, mnl->csg_index_array,
		 mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);

    if( mnl->solver_params.external_inverter == NO_EXT_INV && mnl->solver == BICGSTAB ){
	    fprintf(stderr, "Bicgstab two-step solve not implemented using tmLQCD-native solvers, using CG instead! (det_monomial.c)\n");
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->pf, mnl->solver_params, mnl->maxiter, mnl->forceprec, 
                                     g_relative_precision_flag, VOLUME/2, mnl->Qsq, CG);
    } else {
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->pf, mnl->solver_params, mnl->maxiter, mnl->forceprec, 
                                     g_relative_precision_flag, VOLUME/2, mnl->Qsq, mnl->solver);
    }

    chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    
    /* Y_o -> w_fields[0]  */
    tm_stopwatch_push(&g_timers);
    mnl->Qm(mnl->w_fields[0], mnl->w_fields[1]);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Qm");
    
    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X_e */
    tm_stopwatch_push(&g_timers);
    H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1.);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "H_eo_tm_inv_psi");
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor); 
    
    /* to get the even sites of Y_e */
    tm_stopwatch_push(&g_timers);
    H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "H_eo_tm_inv_psi");
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf, mnl->forcefactor);
  } 
  else {
    /*********************************************************************
     * non even/odd version
     * 
     * This term is det(Q^2 + \mu_1^2)
     *
     *********************************************************************/
    
    if((mnl->solver == CG) || (mnl->solver == MIXEDCG) || (mnl->solver == RGMIXEDCG) || (mnl->solver == MG)) {
      /* Invert Q_{+} Q_{-} */
      /* X -> w_fields[1] */
      chrono_guess(mnl->w_fields[1], mnl->pf, mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_pm_psi);
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->pf, mnl->solver_params, 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, &Q_pm_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME);

      /* Y -> w_fields[0]  */
      tm_stopwatch_push(&g_timers);
      Q_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Q_minus_psi");
      
    }
    else {
      /* Invert first Q_+ */
      /* Y -> w_fields[0]  */
      chrono_guess(mnl->w_fields[0], mnl->pf, mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_plus_psi);
      mnl->iter1 += solve_degenerate(mnl->w_fields[0], mnl->pf, mnl->solver_params, 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, &Q_plus_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[0], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      
      /* Now Q_- */
      /* X -> w_fields[1] */
      
      chrono_guess(mnl->w_fields[1], mnl->w_fields[0], mnl->csg_field2, 
		   mnl->csg_index_array2, mnl->csg_N2, mnl->csg_n2, VOLUME/2, &Q_minus_psi);
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->w_fields[0], mnl->solver_params, 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, &Q_minus_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
        
    }
    
    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(mnl->w_fields[0], mnl->w_fields[1], hf, mnl->forcefactor);
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return;
}


void det_heatbath(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  
  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);
  g_mu = mnl->mu;
  g_kappa = mnl->kappa;
  boundary(g_kappa);
  
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;

  if(mnl->even_odd_flag) {
    tm_stopwatch_push(&g_timers);
    random_spinor_field_eo(mnl->w_fields[0], mnl->rngrepro, RN_GAUSS);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME/2, 1);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "random_energy0");

    tm_stopwatch_push(&g_timers);
    mnl->Qp(mnl->pf, mnl->w_fields[0]);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Qp");

    
    chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    if(mnl->solver != CG) {
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2, 
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
    }
  }
  else {
    tm_stopwatch_push(&g_timers);
    random_spinor_field_lexic(mnl->w_fields[0], mnl->rngrepro,RN_GAUSS);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME, 1);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "random_energy0");

    tm_stopwatch_push(&g_timers);
    Q_plus_psi(mnl->pf, mnl->w_fields[0]);
    tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Q_plus_psi");
    
    chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME);
    if(mnl->solver != CG) {
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2, 
			  mnl->csg_N2, &mnl->csg_n2, VOLUME);
    }
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called det_heatbath for id %d energey %f\n", id, mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return;
}


double det_acc(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  int save_sloppy = g_sloppy_precision_flag;
  
  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);
  g_mu = mnl->mu;
  g_kappa = mnl->kappa;
  boundary(mnl->kappa);
  
  if(mnl->even_odd_flag) {
    g_sloppy_precision_flag = 0;
    if( mnl->solver == MG ){
      chrono_guess(mnl->w_fields[0], mnl->pf, mnl->csg_field, mnl->csg_index_array,
      	     mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qp);
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->pf, mnl->solver_params, mnl->maxiter,
      			      mnl->accprec, g_relative_precision_flag, VOLUME/2, mnl->Qp, mnl->solver);
      tm_stopwatch_push(&g_timers);
      /* Compute the energy contr. from second field */
      mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME/2, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "energy1");
    } else {
      chrono_guess(mnl->w_fields[0], mnl->pf, mnl->csg_field, mnl->csg_index_array,
      mnl->csg_N, mnl->csg_n, VOLUME/2, mnl->Qsq);
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->pf, mnl->solver_params, mnl->maxiter, 
      mnl->accprec, g_relative_precision_flag,VOLUME/2, mnl->Qsq, mnl->solver);
      tm_stopwatch_push(&g_timers);
      mnl->Qm(mnl->w_fields[1], mnl->w_fields[0]);
      /* Compute the energy contr. from first field */
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Qm");

      tm_stopwatch_push(&g_timers);
      mnl->energy1 = square_norm(mnl->w_fields[1], VOLUME/2, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "energy1");
    }
    g_sloppy_precision_flag = save_sloppy;
  }
  else {
    if((mnl->solver == CG) || (mnl->solver == MIXEDCG) || (mnl->solver == RGMIXEDCG)) {
      chrono_guess(mnl->w_fields[1], mnl->pf, mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_pm_psi);
      mnl->iter0 += solve_degenerate(mnl->w_fields[1], mnl->pf, mnl->solver_params, mnl->maxiter, 
                                     mnl->accprec, g_relative_precision_flag, 
			  VOLUME, &Q_pm_psi, mnl->solver);
      
      tm_stopwatch_push(&g_timers);
      Q_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
      /* Compute the energy contr. from first field */
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "Q_minus_psi");
      
      tm_stopwatch_push(&g_timers);
      mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "energy1");
    }
    else {
      chrono_guess(mnl->w_fields[0], mnl->pf, mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_plus_psi);
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->pf, mnl->solver_params, 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME,  &Q_plus_psi, mnl->solver);
      tm_stopwatch_push(&g_timers);
      mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, current_mnl, "energy1");
    }
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called det_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy1 - mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return(mnl->energy1 - mnl->energy0);
}
