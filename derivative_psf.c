/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "linsolve.h"
#include "deriv_Sb.h"
#include "gamma.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "Hopping_Matrix.h"
#include "solver/chrono_guess.h"
#include "derivative_psf.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

void derivative_psf(const int nr, const int set_zero) {
  int i,mu;
#ifdef _KOJAK_INST
#pragma pomp inst begin(derivativepsf)
#endif
  if(set_zero == 1) {
    for(i=0;i<(VOLUME+RAND);i++){ 
      for(mu=0;mu<4;mu++){ 
	_zero_su3adj(df0[i][mu]);
      }
    }
  }

  if(nr == 0) {
    /*********************************************************************
     * 
     * This term is det(Q^2 + \mu_1^2)
     * g_mu1 is set according to the number of psf in use
     *
     *********************************************************************/

    g_mu = g_mu1;
    if(ITER_MAX_BCG == 0 || (fabs(g_mu) > 0) || (g_nr_of_psf != nr+1)) {
      /* If CG is used anyhow */
      /*       gamma5(spionr_field[DUM_DERI+1], g_spinor_field[first_psf], VOLUME/2); */
      
      /* Invert Q_{+} Q_{-} */
      /* X_o -> DUM_DERI+1 */
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[first_psf], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      count00 += solve_cg(DUM_DERI+1, first_psf, g_eps_sq_force1, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
      /*       assign(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+4], VOLUME/2); */
      /* Y_o -> DUM_DERI  */
      Qtm_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
    }
    else{
      /*contributions from field 0 -> first_psf*/
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[first_psf], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      count00 += bicg(DUM_DERI, first_psf, g_eps_sq_force1, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);

      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		   g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
      count01 += bicg(DUM_DERI+1, DUM_DERI, g_eps_sq_force1, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1],
			  g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME/2);
      g_mu = -g_mu;   
    }
  }
  else if(nr == 1) {
    /*********************************************************************
     * 
     * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
     * g_mu1 and g_mu2 are set according to the number of psf in use
     *
     *********************************************************************/
    /* First term coming from the second field */
    /* Multiply with W_+ */
    g_mu = g_mu1;	
    Qtm_plus_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[second_psf]);
    g_mu = g_mu2;
    if(ITER_MAX_BCG == 0 || (fabs(g_mu) > 0) || (g_nr_of_psf != nr+1)) {
      /* If CG is used anyhow */
      /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], VOLUME/2); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      count10 += solve_cg(DUM_DERI+1, DUM_DERI+2, g_eps_sq_force2, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
      /* Y_W -> DUM_DERI  */
      Qtm_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
    }
    else{
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      count10 += bicg(DUM_DERI, DUM_DERI+2, g_eps_sq_force2, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);


      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		   g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
      count11 += bicg(DUM_DERI+1,DUM_DERI, g_eps_sq_force2, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1],
			  g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME/2);
      g_mu = -g_mu;   
    }

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
    g_mu = g_mu1;
    
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(g_spinor_field[DUM_DERI], -1., g_spinor_field[second_psf], VOLUME/2);
    g_mu = g_mu1;
  }

  if(nr == 2) {
    /*********************************************************************
     * 
     * This term is det((Q^2 + \mu_2^2)/(Q^2 + \mu_3^2))
     * g_mu2 and g_mu3 are set according to the number of psf in use
     *
     *********************************************************************/
    /* First term coming from the third field */
    /* Multiply with W_+ */
    g_mu = g_mu2;	
    Qtm_plus_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[third_psf]);
    g_mu = g_mu3;
    if(ITER_MAX_BCG == 0 || (fabs(g_mu) > 0) || (g_nr_of_psf != nr+1)) {
      /* If CG is used anyhow */
      /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2]); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      count20 += solve_cg(DUM_DERI+1, DUM_DERI+2, g_eps_sq_force3, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
      /* Y_W -> DUM_DERI  */
      Qtm_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
    }
    else {
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr],
		   g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
      count20 += bicg(DUM_DERI, DUM_DERI+2, g_eps_sq_force3, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			  g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);

      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		   g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
      gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
      count21 += bicg(DUM_DERI+1,DUM_DERI, g_eps_sq_force3, g_relative_precision_flag);
      chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1],
			  g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME/2);
      g_mu = -g_mu;   
    }

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
    g_mu = g_mu1;

    /* Second term coming from the third field */
    /* The sign is opposite!! */
    mul_r( g_spinor_field[DUM_DERI], -1., g_spinor_field[third_psf], VOLUME/2);
    g_mu = g_mu2;
  }

  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X */
  H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1], EO, -1.);
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
  
  /* to get the even sites of Y */
  H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
  g_mu = g_mu1;
#ifdef _KOJAK_INST
#pragma pomp inst end(derivativepsf)
#endif
}

