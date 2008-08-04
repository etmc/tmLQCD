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
#include "solver/solver.h"
#include "deriv_Sb.h"
#include "tm_operators.h"
#include "chebyshev_polynomial.h"
#include "Nondegenerate_Matrix.h"
#include "Hopping_Matrix.h"
#include "phmc.h"
#include "Nondegenerate_Matrix.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "reweighting_factor_nd.h"
#include "monomial.h"
#include "ndpoly_monomial.h"

extern int phmc_exact_poly;

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void ndpoly_derivative(const int id) {
  int j, k;
  monomial * mnl = &monomial_list[id];

  (*mnl).forcefactor = -phmc_Cpol*phmc_invmaxev;

  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */

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
  } 
  else if(g_epsbar == 0.0) {
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    for(k = 1; k < (phmc_dop_n_cheby-1); k++) {
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[k],
			    g_chi_up_spinor_field[k-1], 
			    phmc_root[k-1]);
    }
    assign(g_chi_up_spinor_field[phmc_dop_n_cheby],
	   g_chi_up_spinor_field[phmc_dop_n_cheby-2], VOLUME/2);

    for(j = (phmc_dop_n_cheby-1); j >= 1; j--) {
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
    Normalisation by the largest  EW  is done in update_momenta
    using mnl->forcefactor
  */ 
}


void ndpoly_heatbath(const int id) {
  int j;
  double temp;
  monomial * mnl = &monomial_list[id];

  printf("called ndpoly_heatbath with g_running_phmc = %d\n", g_running_phmc);
  (*mnl).energy0 = 0.;
  random_spinor_field(g_chi_up_spinor_field[0], VOLUME/2, (*mnl).rngrepro);
  (*mnl).energy0 = square_norm(g_chi_up_spinor_field[0], VOLUME/2);

  if(g_epsbar!=0.0 || phmc_exact_poly == 0){
    random_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2, (*mnl).rngrepro);
     (*mnl).energy0 += square_norm(g_chi_dn_spinor_field[0], VOLUME/2);
  } 
  else {
     zero_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2);
  }

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Here comes the computation of H_old with \n \n");
    printf("PHMC: First: random spinors and their norm  \n ");
    printf("PHMC: OLD Ennergy UP %e \n", (*mnl).energy0);
    printf("PHMC: OLD Energy  DN + UP %e \n\n", (*mnl).energy0);
  }

  if(phmc_exact_poly==0){
    QNon_degenerate(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
		    g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0]);
 
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], 
			phmc_root[phmc_dop_n_cheby-2+j]);
    }
    Poly_tilde_ND(g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], phmc_ptilde_cheby_coef, 
		  phmc_ptilde_n_cheby, g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1]);
  } 
  else if( phmc_exact_poly==1 && g_epsbar!=0.0) {
    /* Attention this is Q * tau1, up/dn are exchanged in the input spinor  */
    /* this is used as an preconditioner */
    QNon_degenerate(g_chi_up_spinor_field[1],g_chi_dn_spinor_field[1],
		    g_chi_dn_spinor_field[0],g_chi_up_spinor_field[0]);

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

    /* solve Q*tau1*P(Q^2) *x=y */
    cg_her_nd(g_chi_up_spinor_field[1],g_chi_dn_spinor_field[1],
	      g_chi_up_spinor_field[0],g_chi_dn_spinor_field[0],
	      1000,1.e-16,0,VOLUME/2, Qtau1_P_ND  ,0,1);

    /*  phi= Bdagger phi  */
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1],
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0],
			phmc_root[phmc_dop_n_cheby-2+j]);
    }

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
  } 
  else if(phmc_exact_poly==1 && g_epsbar==0.0) {
    Qtm_pm_psi(g_chi_up_spinor_field[1], g_chi_up_spinor_field[0]);

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);

    /* solve (Q+)*(Q-)*P((Q+)*(Q-)) *x=y */
    cg_her(g_chi_up_spinor_field[1], g_chi_up_spinor_field[0],
             1000,1.e-16,0,VOLUME/2, Qtm_pm_Ptm_pm_psi  ,0,1);

    /*  phi= Bdagger phi  */
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    phmc_root[phmc_dop_n_cheby-2+j]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
  }

  assign(g_chi_up_copy, g_chi_up_spinor_field[0], VOLUME/2);
  assign(g_chi_dn_copy, g_chi_dn_spinor_field[0], VOLUME/2);

  temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2);
  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Then: evaluate Norm of pseudofermion heatbath BHB \n ");
    printf("PHMC: Norm of BHB up squared %e \n", temp);
  }

  if(g_epsbar!=0.0 || phmc_exact_poly==0) 
    temp += square_norm(g_chi_dn_spinor_field[0], VOLUME/2);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)){
    printf("PHMC: Norm of BHB up + BHB dn squared %e \n\n", temp);
  }
  return;
}


double ndpoly_acc(const int id) {
  int j, ij=0;
  double temp, sgn, fact, Diff;
  double Ener[8];
  double factor[8];
  monomial * mnl = &monomial_list[id];
  
  (*mnl).energy1 = 0.;
  Ener[0] = 0;
  factor[0] = 1.0;
  for(j = 1; j < 8; j++){
    factor[j] = j*factor[j-1];
    Ener[j] = 0;
  }
  /* IF PHMC */

  /* This is needed if we consider only "1" in eq. 9 */
  assign(g_chi_up_spinor_field[1], g_chi_up_copy, VOLUME/2);
  assign(g_chi_dn_spinor_field[1], g_chi_dn_copy,  VOLUME/2);

  if(phmc_exact_poly==0) {
    for(j = 1; j <= (phmc_dop_n_cheby-1); j++) {
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
      
      /* Change this name !!*/
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
			   g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], 
			   phmc_root[j-1]);
    }
  
    ij=0;
    assign(g_chi_up_spinor_field[ij], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[ij], g_chi_dn_spinor_field[1], VOLUME/2);

    temp = square_norm(g_chi_up_spinor_field[ij], VOLUME/2);
    Ener[ij] = temp;

    temp = square_norm(g_chi_dn_spinor_field[ij], VOLUME/2);
    Ener[ij] += temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: Here comes the computation of H_new with \n \n");

      printf("PHMC: At j=%d  P+HMC Final Energy %e \n", ij, (*mnl).energy1+Ener[ij]);
      printf("PHMC: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[ij]);
    }
    
    /* Here comes the loop for the evaluation of A, A^2, ...  */
    for(j = 1; j < 1; j++){ /* To omit corrections just set  j<1 */
      
      if(j % 2){ /*  Chi[j] = ( Qdag P  Ptilde ) Chi[j-1]  */ 
	Poly_tilde_ND(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		      phmc_ptilde_cheby_coef, phmc_ptilde_n_cheby, 
		      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
	QdaggerQ_poly(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		      phmc_dop_cheby_coef, phmc_dop_n_cheby, 
		      g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j]);
	QdaggerNon_degenerate(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
      }
      else { /*  Chi[j] = ( Ptilde P Q ) Chi[j-1]  */ 
	QNon_degenerate(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
	QdaggerQ_poly(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		      phmc_dop_cheby_coef, phmc_dop_n_cheby, g_chi_up_spinor_field[j], 
		      g_chi_dn_spinor_field[j]);
	Poly_tilde_ND(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		      phmc_ptilde_cheby_coef, phmc_ptilde_n_cheby, 
		      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
      }

      Ener[j] = Ener[j-1] + Ener[0];
      sgn = -1.0;
      for(ij = 1; ij < j; ij++){
	fact = factor[j] / (factor[ij] * factor[j-ij]);
	if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
	  printf("PHMC: Here  j=%d  and  ij=%d   sign=%f  fact=%f \n", j ,ij, sgn, fact);
	}
	Ener[j] += sgn*fact*Ener[ij];
	sgn = -sgn;
      }
      temp = square_norm(g_chi_up_spinor_field[j], VOLUME/2);
      temp += square_norm(g_chi_dn_spinor_field[j], VOLUME/2);
      if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
	printf("PHMC: Here  j=%d   sign=%f  temp=%e \n", j, sgn, temp);
      }

      Ener[j] += sgn*temp;

      Diff = fabs(Ener[j] - Ener[j-1]);
      if((g_proc_id == g_stdio_proc) && (g_debug_level > 0)) {
	printf("PHMC: Correction aftern %d steps: %e \n", j, Diff);
      }

      if(Diff < g_acc_Hfin) {
	if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
	  printf("PHMC: At j = %d  PHMC Only Final Energy %e \n", j, Ener[j]);
	}
	break;
      }
    }
    (*mnl).energy1 += Ener[ij];  /* this is quite sticky */
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[ij], (*mnl).energy1);
    }
  } 
  else if(phmc_exact_poly==1 && g_epsbar!=0.0) {
    /* B(Q*tau1) */
    for(j = 1; j <= (phmc_dop_n_cheby-1); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1],
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0],
			phmc_root[j-1]);
    }

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

    temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2);
    Ener[0] = temp;

    temp = square_norm(g_chi_dn_spinor_field[0], VOLUME/2);
    Ener[0] += temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      ij=0;
      printf("PHMC: Here comes the computation of H_new with \n \n");
      printf("PHMC: At j=%d  P+HMC Final Energy %e \n", ij, (*mnl).energy1+Ener[0]);
      printf("PHMC: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    (*mnl).energy1 += Ener[0];
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[0], (*mnl).energy1);
    }
  } 
  else if(phmc_exact_poly == 1 && g_epsbar == 0.0) {
    for(j = 1; j < (phmc_dop_n_cheby); j++) {
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    phmc_root[j-1]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);

    temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2);
    Ener[0] = temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: Here comes the computation of H_new with \n \n");
      printf("PHMC: At j=%d  P+HMC Final Energy %e \n", ij, (*mnl).energy1+Ener[0]);
      printf("PHMC: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    (*mnl).energy1 += Ener[0];
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[0], (*mnl).energy1);
    }
  }
  /* END IF PHMC */
  return(mnl->energy1 - mnl->energy0);
}
