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
#include "su3spinor.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
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
#include "derivative_psf.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

void derivative_psf(const int nr, const int set_zero) 
{
  int i, mu, x;
  extern su3 ** g_stout_force_field;
  extern su3 ** g_gauge_field_saved;

#ifdef _KOJAK_INST
#pragma pomp inst begin(derivativepsf)
#endif
  if(use_stout_flag == 1)
  {

    /*
     *  save unsmeared gauge field
     */
    for(x = 0; x < VOLUME; x++) 
      for(mu = 0; mu < 4; mu++)
      {
        _su3_assign(g_gauge_field_saved[x][mu], g_gauge_field[x][mu]);
      }
    stout_smear_gauge_field(stout_rho , stout_no_iter);
  }

  if(set_zero == 1) 
  {
    for(i=0;i<(VOLUME+RAND);i++)
    { 
      for(mu=0;mu<4;mu++)
      { 
        _zero_su3adj(df0[i][mu]);
      }
    }
  }

  if(even_odd_flag)
  {
    if(nr == 0) 
    {
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
        count00 += solve_cg(g_spinor_field[DUM_DERI+1], g_spinor_field[first_psf], 
			    g_eps_sq_force1, g_relative_precision_flag);
        chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr],
            g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
        /*       assign(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+4], VOLUME/2); */
        /* Y_o -> DUM_DERI  */
        Qtm_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
      }
      else
      {
        /*contributions from field 0 -> first_psf*/
        /* Invert first Q_+ */
        /* Y_o -> DUM_DERI  */
        chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[first_psf], g_csg_field[nr], g_csg_index_array[nr],
		     g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Qtm_pm_psi);
        gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
        count00 += bicg(g_spinor_field[DUM_DERI], g_spinor_field[first_psf]
			, g_eps_sq_force1, g_relative_precision_flag);
        chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			    g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);

        /* Now Q_- */
        /* X_o -> DUM_DERI+1 */
        g_mu = -g_mu;
        chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		     g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
        gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
        count01 += bicg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], 
			g_eps_sq_force1, g_relative_precision_flag);
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
	count10 += solve_cg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], 
			    g_eps_sq_force2, g_relative_precision_flag);
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
	count10 += bicg(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], 
			g_eps_sq_force2, g_relative_precision_flag);
	chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			    g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
	
	
	/* Now Q_- */
	/* X_o -> DUM_DERI+1 */
	g_mu = -g_mu;
	chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		     g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
	gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
	count11 += bicg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], 
			g_eps_sq_force2, g_relative_precision_flag);
	chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1],
			    g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME/2);
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
        count20 += solve_cg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], 
			    g_eps_sq_force3, g_relative_precision_flag);
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
        count20 += bicg(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], 
			g_eps_sq_force3, g_relative_precision_flag);
        chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr],
			    g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME/2);
	
        /* Now Q_- */
        /* X_o -> DUM_DERI+1 */
        g_mu = -g_mu;
        chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1],
		     g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME/2, &Qtm_pm_psi);
        gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
        count21 += bicg(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], 
			g_eps_sq_force3, g_relative_precision_flag);
        chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1],
			    g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME/2);
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
    deriv_Sb(OE, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2]); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+1]);
  } 
  else {
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
        /* X -> DUM_DERI+1 */
        count00 += cg_her(g_spinor_field[DUM_DERI+1], g_spinor_field[first_psf], ITER_MAX_CG, g_eps_sq_force1, g_relative_precision_flag, VOLUME, &Q_pm_psi, 0, 1);
        /*       assign(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+4], VOLUME/2); */
        /* Y -> DUM_DERI  */

        Q_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
        /*Q_minus_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+1]);*/

      }
      else {
        /*contributions from field 0 -> first_psf*/
        /* Invert first Q_+ */
        /* Y_o -> DUM_DERI  */

        /*chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[first_psf], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME, &Q_pm_psi);
          gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);*/
        count00 += bicgstab_complex(g_spinor_field[DUM_DERI], g_spinor_field[first_psf], 1000, g_eps_sq_force1, g_relative_precision_flag, VOLUME,  Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME);*/

        /* Now Q_- */
        /* X_o -> DUM_DERI+1 */
        g_mu = -g_mu;
        /*chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME, &Q_pm_psi);
        gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME);*/
        count01 += bicgstab_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], 1000, g_eps_sq_force1, g_relative_precision_flag, VOLUME, Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME);*/
        g_mu = -g_mu;   
      }
    }
    if(nr == 1) {
      /*********************************************************************
       * 
       * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
       * g_mu1 and g_mu2 are set according to the number of psf in use
       *
       *********************************************************************/
      /* First term coming from the second field */
      /* Multiply with W_+ */

      g_mu = g_mu1;	
      Q_plus_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[second_psf]);
      g_mu = g_mu2;
      if(ITER_MAX_BCG == 0 || (fabs(g_mu) > 0) || (g_nr_of_psf != nr+1)) 
      {
        /* If CG is used anyhow */
        /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], VOLUME/2); */
        /* Invert Q_{+} Q_{-} */
        /* X_W -> DUM_DERI+1 */

        /*chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr],
          g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Q_pm_psi);*/
        count10 += cg_her(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], 
			  ITER_MAX_CG, g_eps_sq_force2, g_relative_precision_flag, 
			  VOLUME, &Q_pm_psi, 0, 1);
        /*chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME);*/
        /* Y_W -> DUM_DERI  */
        Q_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
      }
      else {
        /* Invert first Q_+ */
        /* Y_o -> DUM_DERI  */
        /*chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME, &Qtm_pm_psi);
          gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);*/
        count10 += bicgstab_complex(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], 
				    1000, g_eps_sq_force2, g_relative_precision_flag, 
				    VOLUME, Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME);*/


        /* Now Q_- */
        /* X_o -> DUM_DERI+1 */
        g_mu = -g_mu;
        /*chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME, &Qtm_pm_psi);*/
	gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME);
	count11 += bicgstab_complex(g_spinor_field[DUM_DERI+1],g_spinor_field[DUM_DERI], 
				    1000, g_eps_sq_force2, g_relative_precision_flag, 
				    VOLUME, Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME);*/
        g_mu = -g_mu;   
      }
      
      /* \delta Q sandwitched by Y^\dagger and X */
      deriv_Sb_D_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]); 
      
      g_mu = g_mu1;
      
      /* Second term coming from the second field */
      /* The sign is opposite!! */
      mul_r(g_spinor_field[DUM_DERI], -1., g_spinor_field[second_psf], VOLUME);
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
      Q_plus_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[third_psf]);
      g_mu = g_mu3;
      if(ITER_MAX_BCG == 0 || (fabs(g_mu) > 0) || (g_nr_of_psf != nr+1)) {
        /* If CG is used anyhow */
        /*       gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2]); */
        /* Invert Q_{+} Q_{-} */
        /* X_W -> DUM_DERI+1 */
        /*chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME/2, &Q_pm_psi);*/
        count20 += cg_her(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], 
			  1000, g_eps_sq_force3, g_relative_precision_flag, 
			  VOLUME, &Q_pm_psi, 0, 1);
        /*chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME);*/
        /* Y_W -> DUM_DERI  */
        Q_minus_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
      }
      else {
        /* Invert first Q_+ */
        /* Y_o -> DUM_DERI  */
        /*chrono_guess(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], g_csg_N[2*nr+1], VOLUME, &Q_pm_psi);*/
        gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME);
        count20 += bicgstab_complex(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+2], 
				    1000, g_eps_sq_force3, g_relative_precision_flag, 
				    VOLUME, Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI], g_csg_field[nr], g_csg_index_array[nr], g_csg_N[2*nr], &g_csg_N[2*nr+1], VOLUME);*/
	
        /* Now Q_- */
        /* X_o -> DUM_DERI+1 */
        g_mu = -g_mu;
        /*chrono_guess(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], g_csg_N[2*nr+3], VOLUME, &Q_pm_psi);*/
        gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME);
        count21 += bicgstab_complex(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI], 
				    1000, g_eps_sq_force3, g_relative_precision_flag, 
				    VOLUME, Q_minus_psi);
        /*chrono_add_solution(g_spinor_field[DUM_DERI+1], g_csg_field[nr+1], g_csg_index_array[nr+1], g_csg_N[2*nr+2], &g_csg_N[2*nr+3], VOLUME);*/
        g_mu = -g_mu;   
      }

      deriv_Sb_D_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
      g_mu = g_mu1;

      /* Second term coming from the third field */
      /* The sign is opposite!! */
      mul_r( g_spinor_field[DUM_DERI], -1., g_spinor_field[third_psf], VOLUME);
      g_mu = g_mu2;

    }

    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  }
  g_mu = g_mu1;

  if(use_stout_flag == 1)
  {
    /*
     *  now we iterate the force field (\Sigma in hep-lat/0311018) 
     *  according to eqtn(75) in hep-lat/0311018
     *  for this we need the force terms as explicit matrices
     */
    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < 4; mu++)
      {
        _make_su3(g_stout_force_field[x][mu], df0[x][mu]);
      }
    stout_smear_force();

    for(x = 0; x < VOLUME; x++)
      for(mu = 0; mu < 4; mu++)
      {
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

    printf("df0 after = %f %f %f %f %f %f %f %f\n", df0[0][0].d1, df0[0][0].d2, df0[0][0].d3, df0[0][0].d4, df0[0][0].d5, df0[0][0].d6, df0[0][0].d7, df0[0][0].d8);

    /*
     *  restore unsmeared gauge field
     */
    for(x = 0; x < VOLUME; x++) 
      for(mu = 0; mu < 4; mu++)
      {
        _su3_assign(g_gauge_field[x][mu], g_gauge_field_saved[x][mu]);
      }
  }

#ifdef _KOJAK_INST
#pragma pomp inst end(derivativepsf)
#endif
}
