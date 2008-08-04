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
#include "linalg_eo.h"
#include "start.h"
#include "linsolve.h"
#include "deriv_Sb.h"
#include "tm_operators.h"
#include "chebyshev_polynomial.h"
#include "Nondegenerate_Matrix.h"
#include "Hopping_Matrix.h"
#include "phmc.h"
#include "derivative_nondegenerate.h"


extern int phmc_exact_poly;

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void derivative_nondegenerate() {

  int i, mu, j, k;

  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */
  /* Re-initialize df0 */
  for(i=0;i<(VOLUME+RAND);i++) { 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }

  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    for(k = 1; k < (phmc_dop_n_cheby-1); k++) {
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[k], g_chi_dn_spinor_field[k], 
			   g_chi_up_spinor_field[k-1], g_chi_dn_spinor_field[k-1], 
			   phmc_root[k-1]);
    }
    
    
    /* Here comes the remaining fields  chi_k ; k=n,...,2n-1  */
    /*They are evaluated step-by-step overwriting the same field (phmc_dop_n_cheby)*/
    
    assign(g_chi_up_spinor_field[phmc_dop_n_cheby], g_chi_up_spinor_field[phmc_dop_n_cheby-2], VOLUME/2);
    assign(g_chi_dn_spinor_field[phmc_dop_n_cheby], g_chi_dn_spinor_field[phmc_dop_n_cheby-2], VOLUME/2);
    
    for(j=(phmc_dop_n_cheby-1); j>=1; j--) {
      
      assign(g_chi_up_spinor_field[phmc_dop_n_cheby-1], g_chi_up_spinor_field[phmc_dop_n_cheby], VOLUME/2);
      assign(g_chi_dn_spinor_field[phmc_dop_n_cheby-1], g_chi_dn_spinor_field[phmc_dop_n_cheby], VOLUME/2);
      
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[phmc_dop_n_cheby], g_chi_dn_spinor_field[phmc_dop_n_cheby], 
			   g_chi_up_spinor_field[phmc_dop_n_cheby-1], g_chi_dn_spinor_field[phmc_dop_n_cheby-1], 
			   phmc_root[2*phmc_dop_n_cheby-3-j]);
      
      /*     assign(g_spinor_field[DUM_DERI+4], g_chi_up_spinor_field[phmc_dop_n_cheby], VOLUME/2); */
      /*     assign(g_spinor_field[DUM_DERI+5], g_chi_dn_spinor_field[phmc_dop_n_cheby], VOLUME/2); */
      
      /*     assign(g_spinor_field[DUM_DERI+2], g_chi_up_spinor_field[j-1], VOLUME/2); */
      /*     assign(g_spinor_field[DUM_DERI+3], g_chi_dn_spinor_field[j-1], VOLUME/2); */
      
      /* Get the even parts of the  (j-1)th  chi_spinors */
      H_eo_ND(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], 
	      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], EO);
      
      /* \delta M_eo sandwitched by  chi[j-1]_e^\dagger  and  chi[2N-j]_o */
      deriv_Sb(EO, g_spinor_field[DUM_DERI], g_chi_up_spinor_field[phmc_dop_n_cheby]);      /* UP */
      deriv_Sb(EO, g_spinor_field[DUM_DERI+1], g_chi_dn_spinor_field[phmc_dop_n_cheby]);    /* DN */
      
      /* Get the even parts of the  (2N-j)-th  chi_spinors */
      H_eo_ND(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], 
	      g_chi_up_spinor_field[phmc_dop_n_cheby], g_chi_dn_spinor_field[phmc_dop_n_cheby], EO);
      
      /* \delta M_oe sandwitched by  chi[j-1]_o^\dagger  and  chi[2N-j]_e */
      deriv_Sb(OE, g_chi_up_spinor_field[j-1], g_spinor_field[DUM_DERI]);
      deriv_Sb(OE, g_chi_dn_spinor_field[j-1], g_spinor_field[DUM_DERI+1]);
    }
  } else  if(g_epsbar==0.0){

    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */

    for(k = 1; k < (phmc_dop_n_cheby-1); k++) {
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[k],
			    g_chi_up_spinor_field[k-1], 
			    phmc_root[k-1]);
    }


    assign(g_chi_up_spinor_field[phmc_dop_n_cheby],
	   g_chi_up_spinor_field[phmc_dop_n_cheby-2], VOLUME/2);

    for(j=(phmc_dop_n_cheby-1); j>=1; j--) {

      assign(g_chi_up_spinor_field[phmc_dop_n_cheby-1],
	     g_chi_up_spinor_field[phmc_dop_n_cheby], VOLUME/2);

      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[phmc_dop_n_cheby], 
			   g_chi_up_spinor_field[phmc_dop_n_cheby-1],
			   phmc_root[2*phmc_dop_n_cheby-3-j]);

      Qtm_minus_psi(g_spinor_field[DUM_DERI+3],g_chi_up_spinor_field[j-1]); 

      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_chi_up_spinor_field[phmc_dop_n_cheby], EO, -1.);
      deriv_Sb(OE, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+2]); 
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], EO, 1.); 
      deriv_Sb(EO, g_spinor_field[DUM_DERI+2], g_chi_up_spinor_field[phmc_dop_n_cheby]);




      Qtm_minus_psi(g_spinor_field[DUM_DERI+3],g_chi_up_spinor_field[phmc_dop_n_cheby]); 


      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2],g_spinor_field[DUM_DERI+3], EO, +1.);
      deriv_Sb(OE, g_chi_up_spinor_field[j-1] , g_spinor_field[DUM_DERI+2]); 
      
      H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_chi_up_spinor_field[j-1], EO, -1.); 
      deriv_Sb(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);

    }

  }

  
  /*
    Normalisation by the largest  EW  is done in update_fermion_momenta
    C.U. this is something worth to be changed...
  */ 
}


