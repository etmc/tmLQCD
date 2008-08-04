/* $Id$ */

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
#include "stout_smear.h"
#include "stout_smear_force.h"
#include "monomial.h"
#include "boundary.h"
#include "detratio_monomial.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

/* think about chronological solver ! */

void detratio_derivative(const int no) {
  int mu, x, saveiter = ITER_MAX_BCG;
  extern su3 ** g_stout_force_field;
  extern su3 ** g_gauge_field_saved;
  monomial * mnl = &monomial_list[no];
  spinor * psf = (*mnl).pf;

  mnl->forcefactor = 1.;

  if(use_stout_flag == 1) {
    /*  save unsmeared gauge field */
    for(x = 0; x < VOLUME; x++) {
      for(mu = 0; mu < 4; mu++) {
        _su3_assign(g_gauge_field_saved[x][mu], g_gauge_field[x][mu]);
      }
    }
    stout_smear_gauge_field(stout_rho , stout_no_iter);
  }

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
    Qtm_plus_psi(g_spinor_field[DUM_DERI+2], psf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(mnl->solver == CG || (fabs(g_mu) > 0)) {
      ITER_MAX_CG = mnl->maxiter;
      /* If CG is used anyhow */
      /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], VOLUME/2); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_pm_psi);
      mnl->iter1 += solve_cg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], mnl->forceprec, 
			     g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      /* Y_W -> DUM_DERI  */
      Qtm_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
    }
    else {
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      ITER_MAX_BCG = mnl->maxiter;
      chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], mnl->csg_field, mnl->csg_index_array,
		   mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      mnl->iter1 += bicg(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], mnl->forceprec, 
			 g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      
      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], mnl->csg_field2, mnl->csg_index_array2,
		   mnl->csg_N2, mnl->csg_n2, VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
      mnl->iter1 += bicg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], mnl->forceprec, 
			 g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
      g_mu = -g_mu;   
    }
    
    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2]); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1]); 
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(g_spinor_field[DUM_DERI], -1., psf, VOLUME/2);

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2]); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1]);
  } 
  else {
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
    Q_plus_psi(g_spinor_field[DUM_DERI+2], psf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(mnl->solver == CG || (fabs(g_mu) > 0)) {/*  || (g_nr_of_psf != nr+1)) { */
      /* If CG is used anyhow */
      /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], VOLUME/2); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      mnl->iter1 += cg_her(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], 
			   mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
			   VOLUME, &Q_pm_psi, 0, 1);
      /* Y_W -> DUM_DERI  */
      Q_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
    }
    else {
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      /*gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);*/
      mnl->iter1 += bicgstab_complex(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, Q_minus_psi);
      
      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME);
      mnl->iter1 += bicgstab_complex(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI], 
				     mnl->maxiter, mnl->forceprec, g_relative_precision_flag, 
				     VOLUME, Q_minus_psi);
      g_mu = -g_mu;   
    }

    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]); 
    
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(g_spinor_field[DUM_DERI], -1., psf, VOLUME);
    
    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  }
  g_mu = g_mu1;
  boundary(g_kappa);

  if(use_stout_flag == 1) {
    /*
     *  now we iterate the force field (\Sigma in hep-lat/0311018) 
     *  according to eqtn(75) in hep-lat/0311018
     *  for this we need the force terms as explicit matrices
     */
    for(x = 0; x < VOLUME; x++) {
      for(mu = 0; mu < 4; mu++) {
        _make_su3(g_stout_force_field[x][mu], df0[x][mu]);
      }
    }
    stout_smear_force();
    
    for(x = 0; x < VOLUME; x++) {
      for(mu = 0; mu < 4; mu++) {
        _trace_lambda(df0[x][mu],g_stout_force_field[x][mu]);
        df0[x][mu].d1 /= -2.0;
        df0[x][mu].d2 /= -2.0;
        df0[x][mu].d3 /= -2.0;
        df0[x][mu].d4 /= -2.0;
        df0[x][mu].d5 /= -2.0;
        df0[x][mu].d6 /= -2.0;
        df0[x][mu].d7 /= -2.0;
        df0[x][mu].d8 /= -2.0;
      }
    }

    printf("df0 after = %f %f %f %f %f %f %f %f\n", df0[0][0].d1, df0[0][0].d2, df0[0][0].d3, df0[0][0].d4, df0[0][0].d5, df0[0][0].d6, df0[0][0].d7, df0[0][0].d8);

    /*
     *  restore unsmeared gauge field
     */
    for(x = 0; x < VOLUME; x++) {
      for(mu = 0; mu < 4; mu++) {
        _su3_assign(g_gauge_field[x][mu], g_gauge_field_saved[x][mu]);
      }
    }
  }
  ITER_MAX_BCG = saveiter;
  return;
}


void detratio_heatbath(const int id) {
  int saveiter = ITER_MAX_BCG;
  monomial * mnl = &monomial_list[id];

  g_mu = mnl->mu;
  boundary(mnl->kappa);
  mnl->csg_n = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;
  if(mnl->even_odd_flag) {
    random_spinor_field(g_spinor_field[4], VOLUME/2, mnl->rngrepro);
    mnl->energy0  = square_norm(g_spinor_field[4], VOLUME/2);

    Qtm_plus_psi(g_spinor_field[3], g_spinor_field[4]);
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    zero_spinor_field(mnl->pf,VOLUME/2);
    if(mnl->solver == CG || fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    ITER_MAX_CG = mnl->maxiter;
    mnl->iter0 += bicg(mnl->pf, g_spinor_field[3], mnl->accprec, g_relative_precision_flag);

    chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array,
			mnl->csg_N, &mnl->csg_n, VOLUME/2);
    if(g_nr_of_psf == 2 && ITER_MAX_BCG > 0 && fabs(g_mu) == 0.) {
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2,
			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2);
    }
  }
  else {
    random_spinor_field(g_spinor_field[4], VOLUME, mnl->rngrepro);
    mnl->energy0 = square_norm(g_spinor_field[4], VOLUME);

    Q_plus_psi(g_spinor_field[3], g_spinor_field[4]);
    g_mu = mnl->mu2;
    boundary(mnl->kappa2);
    zero_spinor_field(mnl->pf,VOLUME);
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    mnl->iter0 += bicgstab_complex(mnl->pf, g_spinor_field[3], mnl->maxiter, mnl->accprec, 
				   g_relative_precision_flag, VOLUME, Q_minus_psi);
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called detratio_heatbath for id %d %d energy %f\n", id, mnl->even_odd_flag, mnl->energy0);
  }
  g_mu = g_mu1;
  boundary(g_kappa);
  ITER_MAX_BCG = saveiter;
  return;
}

double detratio_acc(const int id) {
  monomial * mnl = &monomial_list[id];
  int saveiter = ITER_MAX_BCG;

  g_mu = mnl->mu2;
  boundary(mnl->kappa2);
  if(even_odd_flag) {
    Qtm_plus_psi(g_spinor_field[DUM_DERI+5], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(mnl->solver == CG || fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    ITER_MAX_CG = mnl->maxiter;
    chrono_guess(g_spinor_field[3], g_spinor_field[DUM_DERI+5], mnl->csg_field, mnl->csg_index_array, 
		 mnl->csg_N, mnl->csg_n, VOLUME/2, &Qtm_pm_psi);
    mnl->iter0 += bicg(g_spinor_field[3], g_spinor_field[DUM_DERI+5], mnl->accprec, g_relative_precision_flag); 
    /*     ITER_MAX_BCG = *saveiter_max; */
    /* Compute the energy contr. from second field */
    mnl->energy1 = square_norm(g_spinor_field[3], VOLUME/2);
  }
  else {
    Q_plus_psi(g_spinor_field[DUM_DERI+5], mnl->pf);
    g_mu = mnl->mu;
    boundary(mnl->kappa);
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    mnl->iter0 += bicgstab_complex(g_spinor_field[3], g_spinor_field[DUM_DERI+5], 
				   mnl->maxiter, mnl->accprec, g_relative_precision_flag, VOLUME, Q_minus_psi); 
    /*     ITER_MAX_BCG = *saveiter_max; */
    /* Compute the energy contr. from second field */
    mnl->energy1 = square_norm(g_spinor_field[3], VOLUME);
  }
  g_mu = g_mu1;
  boundary(g_kappa);
  ITER_MAX_BCG = saveiter;
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called detratio_acc for id %d %d energy %f\n", id, mnl->even_odd_flag, mnl->energy1);
  }
  return(mnl->energy1 - mnl->energy0);
}
