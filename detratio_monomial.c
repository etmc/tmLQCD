/***********************************************************************
 *
 * Copyright (C) 1008 Carsten Urbach
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
#include "solver/bicgstab_complex.h"
#include "solver/solver.h"
#include "read_input.h"
#include "smearing/stout.h"
#include "clover_leaf.h"

#include "monomial.h"
#include "boundary.h"
#include "detratio_monomial.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

/* think about chronological solver ! */

void detratio_derivative(const int no, hamiltonian_field_t * const hf) {
  int saveiter = ITER_MAX_BCG;

  monomial * mnl = &monomial_list[no];

  /* This factor 2* a missing factor 2 in trace_lambda */
  mnl->forcefactor = 1.;

  if(mnl->even_odd_flag) {
    /*
     *  this is being run in case there is even/odd preconditioning
     */
    /*********************************************************************
     * 
     * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
     * mu1 and mu2 are set according to the monomial
     *
     *********************************************************************/
    /* First term coming from the second field */
    /* Multiply with W_+ */
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);

    if(mnl->solver != CG) {
      fprintf(stderr, "Bicgstab currently not implemented, using CG instead! (detratio_monomial.c)\n");
    }

    Qtm_plus_psi(mnl->w_fields[2], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    /* Invert Q_{+} Q_{-} */
    /* X_W -> w_fields[1] */
    chrono_guess(mnl->w_fields[1], mnl->w_fields[2], mnl->csg_field, 
		 mnl->csg_index_array, mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_pm_psi);
    mnl->iter1 += cg_her(mnl->w_fields[1], mnl->w_fields[2], mnl->maxiter, 
			 mnl->forceprec, g_relative_precision_flag, VOLUME/2, &Qtm_pm_psi);
    chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    /* Y_W -> w_fields[0]  */
    Qtm_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
    
    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf); 

    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(mnl->w_fields[0], -1., mnl->pf, VOLUME/2);

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf);

  } 
  else { /* no even/odd preconditioning */
    /*********************************************************************
     * 
     * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
     * mu1 and mu2 are set according to the monomial
     *
     *********************************************************************/
    /* First term coming from the second field */
    /* Multiply with W_+ */
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);	
    Q_plus_psi(mnl->w_fields[2], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(mnl->solver == CG) {
      /* If CG is used anyhow */
      /*       gamma5(mnl->w_fields[1], mnl->w_fields[2], VOLUME/2); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> w_fields[1] */
      chrono_guess(mnl->w_fields[1], mnl->w_fields[2], mnl->csg_field, 
		   mnl->csg_index_array, mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_pm_psi);
      mnl->iter1 += cg_her(mnl->w_fields[1], mnl->w_fields[2], 
			   mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
			   VOLUME, &Q_pm_psi);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      
      /* Y_W -> w_fields[0]  */
      Q_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
    }
    else {
      /* Invert first Q_+ */
      /* Y_o -> w_fields[0]  */

      chrono_guess(mnl->w_fields[0], mnl->w_fields[2], mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_plus_psi);
      gamma5(mnl->w_fields[0], mnl->w_fields[0], VOLUME);
      mnl->iter1 += bicgstab_complex(mnl->w_fields[0], mnl->w_fields[2], 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, Q_plus_psi);
      chrono_add_solution(mnl->w_fields[0], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);

      /* Now Q_- */
      /* X_o -> w_fields[1] */
      g_mu = -g_mu;
      chrono_guess(mnl->w_fields[1], mnl->w_fields[0], mnl->csg_field2, 
		   mnl->csg_index_array2, mnl->csg_N2, mnl->csg_n2, VOLUME/2, &Q_minus_psi);
      gamma5(mnl->w_fields[1], mnl->w_fields[1], VOLUME);
      mnl->iter1 += bicgstab_complex(mnl->w_fields[1],mnl->w_fields[0], 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, Q_minus_psi);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
      g_mu = -g_mu;   
    }

    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(mnl->w_fields[0], mnl->w_fields[1], hf); 
    
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(mnl->w_fields[0], -1., mnl->pf, VOLUME);
    
    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(mnl->w_fields[0], mnl->w_fields[1], hf);
  }
  g_mu = g_mu1;
  boundary(g_kappa);

  ITER_MAX_BCG = saveiter;
  return;
}


void detratio_heatbath(const int id, hamiltonian_field_t * const hf) {
  int saveiter = ITER_MAX_BCG;
  monomial * mnl = &monomial_list[id];

  g_mu = mnl->mu;
  boundary(mnl->kappa);
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;
  if(mnl->even_odd_flag) {
    random_spinor_field(mnl->w_fields[0], VOLUME/2, mnl->rngrepro);
    mnl->energy0  = square_norm(mnl->w_fields[0], VOLUME/2, 1);

    Qtm_plus_psi(mnl->w_fields[1], mnl->w_fields[0]);
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    zero_spinor_field(mnl->pf,VOLUME/2);
    if(mnl->solver == CG) ITER_MAX_BCG = 0;
    ITER_MAX_CG = mnl->maxiter;
    mnl->iter0 += bicg(mnl->pf, mnl->w_fields[1], mnl->accprec, g_relative_precision_flag);

    chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    if(mnl->solver != CG) {
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
    }
  }
  else {
    random_spinor_field(mnl->w_fields[0], VOLUME, mnl->rngrepro);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME, 1);

    Q_plus_psi(mnl->w_fields[1], mnl->w_fields[0]);
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    zero_spinor_field(mnl->pf,VOLUME);
    mnl->iter0 += bicgstab_complex(mnl->pf, mnl->w_fields[1], mnl->maxiter, mnl->accprec, 
				   g_relative_precision_flag, VOLUME, Q_plus_psi);
    chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    if(mnl->solver != CG) {
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
    }
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called detratio_heatbath for id %d %d energy %f\n", id, mnl->even_odd_flag, mnl->energy0);
  }
  g_mu = g_mu1;
  boundary(g_kappa);
  ITER_MAX_BCG = saveiter;
  return;
}

double detratio_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  int saveiter = ITER_MAX_BCG;
  int save_sloppy = g_sloppy_precision_flag;

  g_mu = mnl->mu2;
  boundary(mnl->kappa2);
  if(even_odd_flag) {
    Qtm_plus_psi(mnl->w_fields[1], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(mnl->solver == CG) ITER_MAX_BCG = 0;
    ITER_MAX_CG = mnl->maxiter;
    chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array, 
		 mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_plus_psi);
    g_sloppy_precision_flag = 0;    
    mnl->iter0 += bicg(mnl->w_fields[0], mnl->w_fields[1], mnl->accprec, g_relative_precision_flag); 
    g_sloppy_precision_flag = save_sloppy;
    /*     ITER_MAX_BCG = *saveiter_max; */
    /* Compute the energy contr. from second field */
    mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME/2, 1);
  }
  else {
    Q_plus_psi(mnl->w_fields[1], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array, 
		 mnl->csg_N, mnl->csg_n, VOLUME/2, &Q_plus_psi);
    mnl->iter0 += bicgstab_complex(mnl->w_fields[0], mnl->w_fields[1], 
				   mnl->maxiter, mnl->accprec, g_relative_precision_flag, 
				   VOLUME, Q_plus_psi); 
    /*     ITER_MAX_BCG = *saveiter_max; */
    /* Compute the energy contr. from second field */
    mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME, 1);
  }
  g_mu = g_mu1;
  boundary(g_kappa);
  ITER_MAX_BCG = saveiter;
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called detratio_acc for id %d %d dH = %1.4e\n", 
	   id, mnl->even_odd_flag, mnl->energy1 - mnl->energy0);
  }
  return(mnl->energy1 - mnl->energy0);
}
