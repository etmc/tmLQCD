/***********************************************************************
 *
 * Copyright (C) 2008 Thomas Chiarappa, Carsten Urbach
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
#include "hamiltonian_field.h"
#include "boundary.h"
#include "phmc.h"
#include "init_chi_spinor_field.h"
#include "ndpoly_monomial.h"

extern int phmc_exact_poly;
void ndpoly_set_global_parameter(monomial * const mnl);

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void ndpoly_derivative(const int id, hamiltonian_field_t * const hf) {
  int j, k;
  monomial * mnl = &monomial_list[id];


  /* This factor 2 a missing factor 2 in trace_lambda */
  ndpoly_set_global_parameter(mnl);
  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
  }
  else {
    phmc_Cpol = mnl->MDPolyLocNormConst;
  }
  mnl->forcefactor = -2.*phmc_Cpol*phmc_invmaxev;

  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */

  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    /* in  g_chi_up_spinor_field[0] (g_chi_dn_spinor_field[0] we expect */
    /* to find the phi field, the pseudo fermion field                  */
    /* i.e. must be equal to mnl->pf (mnl->pf2)                         */

    assign(g_chi_up_spinor_field[0], mnl->pf, VOLUME/2);
    assign(g_chi_dn_spinor_field[0], mnl->pf2, VOLUME/2);

    for(k = 1; k < (phmc_dop_n_cheby-1); k++) {
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[k], g_chi_dn_spinor_field[k], 
			   g_chi_up_spinor_field[k-1], g_chi_dn_spinor_field[k-1], 
			   mnl->MDPolyRoots[k-1]);
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
			   mnl->MDPolyRoots[2*phmc_dop_n_cheby-3-j]);
      
      /* Get the even parts of the  (j-1)th  chi_spinors */
      H_eo_ND(mnl->w_fields[0], mnl->w_fields[1], 
	      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], EO);
      
      /* \delta M_eo sandwitched by  chi[j-1]_e^\dagger  and  chi[2N-j]_o */
      deriv_Sb(EO, mnl->w_fields[0], g_chi_up_spinor_field[phmc_dop_n_cheby], hf);      /* UP */
      deriv_Sb(EO, mnl->w_fields[1], g_chi_dn_spinor_field[phmc_dop_n_cheby], hf);    /* DN */
      
      /* Get the even parts of the  (2N-j)-th  chi_spinors */
      H_eo_ND(mnl->w_fields[0], mnl->w_fields[1], 
	      g_chi_up_spinor_field[phmc_dop_n_cheby], g_chi_dn_spinor_field[phmc_dop_n_cheby], EO);
      
      /* \delta M_oe sandwitched by  chi[j-1]_o^\dagger  and  chi[2N-j]_e */
      deriv_Sb(OE, g_chi_up_spinor_field[j-1], mnl->w_fields[0], hf);
      deriv_Sb(OE, g_chi_dn_spinor_field[j-1], mnl->w_fields[1], hf);
    }
  } 
  else if(g_epsbar == 0.0) {
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    assign(g_chi_up_spinor_field[0], mnl->pf, VOLUME/2);
    for(k = 1; k < (phmc_dop_n_cheby-1); k++) {
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[k],
			    g_chi_up_spinor_field[k-1], 
			    mnl->MDPolyRoots[k-1]);
    }
    assign(g_chi_up_spinor_field[phmc_dop_n_cheby],
	   g_chi_up_spinor_field[phmc_dop_n_cheby-2], VOLUME/2);

    for(j = (phmc_dop_n_cheby-1); j >= 1; j--) {
      assign(g_chi_up_spinor_field[phmc_dop_n_cheby-1],
	     g_chi_up_spinor_field[phmc_dop_n_cheby], VOLUME/2);

      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[phmc_dop_n_cheby], 
			   g_chi_up_spinor_field[phmc_dop_n_cheby-1],
			   mnl->MDPolyRoots[2*phmc_dop_n_cheby-3-j]);

      Qtm_minus_psi(mnl->w_fields[3],g_chi_up_spinor_field[j-1]); 

      H_eo_tm_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[phmc_dop_n_cheby], EO, -1.);
      deriv_Sb(OE, mnl->w_fields[3], mnl->w_fields[2], hf); 
      
      H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[3], EO, 1.); 
      deriv_Sb(EO, mnl->w_fields[2], g_chi_up_spinor_field[phmc_dop_n_cheby], hf);

      Qtm_minus_psi(mnl->w_fields[3],g_chi_up_spinor_field[phmc_dop_n_cheby]); 

      H_eo_tm_inv_psi(mnl->w_fields[2],mnl->w_fields[3], EO, +1.);
      deriv_Sb(OE, g_chi_up_spinor_field[j-1] , mnl->w_fields[2], hf); 
      
      H_eo_tm_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[j-1], EO, -1.); 
      deriv_Sb(EO, mnl->w_fields[2], mnl->w_fields[3], hf);
    }
  }
  /*
    Normalisation by the largest  EW  is done in update_momenta
    using mnl->forcefactor
  */ 
}


void ndpoly_heatbath(const int id, hamiltonian_field_t * const hf) {
  int j;
  double temp;
  monomial * mnl = &monomial_list[id];

  ndpoly_set_global_parameter(mnl);
  mnl->energy0 = 0.;
  random_spinor_field(g_chi_up_spinor_field[0], VOLUME/2, mnl->rngrepro);
  mnl->energy0 = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);

  if(g_epsbar!=0.0 || phmc_exact_poly == 0) {
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
    random_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2, mnl->rngrepro);
     mnl->energy0 += square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);
  } 
  else {
    phmc_Cpol = mnl->MDPolyLocNormConst;
    zero_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2);
  }

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Here comes the computation of H_old with \n \n");
    printf("PHMC: First: random spinors and their norm  \n ");
    printf("PHMC: OLD Ennergy UP %e \n", mnl->energy0);
    printf("PHMC: OLD Energy  DN + UP %e \n\n", mnl->energy0);
  }

  if(phmc_exact_poly==0){
    QNon_degenerate(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
		    g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0]);
 
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], 
			mnl->MDPolyRoots[phmc_dop_n_cheby-2+j]);
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
	      1000,1.e-16,0,VOLUME/2, Qtau1_P_ND);

    /*  phi= Bdagger phi  */
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
      Q_tau1_min_cconst_ND(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1],
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0],
			mnl->MDPolyRoots[phmc_dop_n_cheby-2+j]);
    }

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
  } 
  else if(phmc_exact_poly==1 && g_epsbar==0.0) {
    Qtm_pm_psi(g_chi_up_spinor_field[1], g_chi_up_spinor_field[0]);

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);

    /* solve (Q+)*(Q-)*P((Q+)*(Q-)) *x=y */
    cg_her(g_chi_up_spinor_field[1], g_chi_up_spinor_field[0],
             1000,1.e-16,0,VOLUME/2, Qtm_pm_Ptm_pm_psi);

    /*  phi= Bdagger phi  */
    for(j = 1; j < (phmc_dop_n_cheby); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    mnl->MDPolyRoots[phmc_dop_n_cheby-2+j]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
  }

  assign(mnl->pf, g_chi_up_spinor_field[0], VOLUME/2);
  assign(mnl->pf2, g_chi_dn_spinor_field[0], VOLUME/2);

  temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);
  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Then: evaluate Norm of pseudofermion heatbath BHB \n ");
    printf("PHMC: Norm of BHB up squared %e \n", temp);
  }

  if(g_epsbar!=0.0 || phmc_exact_poly==0) 
    temp += square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)){
    printf("PHMC: Norm of BHB up + BHB dn squared %e \n\n", temp);
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called ndpoly_heatbath for id %d with g_running_phmc = %d\n", id, g_running_phmc);
  }
  return;
}


double ndpoly_acc(const int id, hamiltonian_field_t * const hf) {
  int j, ij=0;
  double temp, sgn, fact, Diff;
  double Ener[8];
  double factor[8];
  monomial * mnl = &monomial_list[id];
  spinor *up0, *dn0, *up1, *dn1, *dummy;

  ndpoly_set_global_parameter(mnl);
  mnl->energy1 = 0.;
  Ener[0] = 0;
  factor[0] = 1.0;
  for(j = 1; j < 8; j++){
    factor[j] = j*factor[j-1];
    Ener[j] = 0;
  }
  /* IF PHMC */
  up0 = g_chi_up_spinor_field[0];
  up1 = g_chi_up_spinor_field[1];
  dn0 = g_chi_dn_spinor_field[0];
  dn1 = g_chi_dn_spinor_field[1];
  /* This is needed if we consider only "1" in eq. 9 */
  assign(up0, mnl->pf , VOLUME/2);
  assign(dn0, mnl->pf2, VOLUME/2);

  if(phmc_exact_poly==0) {
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
    for(j = 1; j <= (phmc_dop_n_cheby-1); j++) {
      /* Change this name !!*/
      Q_tau1_min_cconst_ND(up1, dn1, up0, dn0, mnl->MDPolyRoots[j-1]);

      dummy = up1; up1 = up0; up0 = dummy;
      dummy = dn1; dn1 = dn0; dn0 = dummy;
      /* result always in up0 and dn0 */
    }
  
    ij=0;
    if(up0 != g_chi_up_spinor_field[ij]) {
      assign(g_chi_up_spinor_field[ij], up0, VOLUME/2);
      assign(g_chi_dn_spinor_field[ij], dn0, VOLUME/2);
    }

    temp = square_norm(g_chi_up_spinor_field[ij], VOLUME/2, 1);
    Ener[ij] = temp;

    temp = square_norm(g_chi_dn_spinor_field[ij], VOLUME/2, 1);
    Ener[ij] += temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: Here comes the computation of H_new with \n \n");

      printf("PHMC: At j=%d  PHMC Final Energy %e \n", ij, mnl->energy1+Ener[ij]);
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
      temp = square_norm(g_chi_up_spinor_field[j], VOLUME/2, 1);
      temp += square_norm(g_chi_dn_spinor_field[j], VOLUME/2, 1);
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
    mnl->energy1 += Ener[ij];  /* this is quite sticky */
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[ij], mnl->energy1);
    }
  } 
  else if(phmc_exact_poly==1 && g_epsbar!=0.0) {
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
    /* B(Q*tau1) */
    for(j = 1; j <= (phmc_dop_n_cheby-1); j++){
      Q_tau1_min_cconst_ND(up1, dn1, up0, dn0, mnl->MDPolyRoots[j-1]);

      dummy = up1; up1 = up0; up0 = dummy;
      dummy = dn1; dn1 = dn0; dn0 = dummy;
      /* result always in up0 and dn0 */
    }
    if(up0 != g_chi_up_spinor_field[0]) {
      assign(g_chi_up_spinor_field[0], up0, VOLUME/2);
      assign(g_chi_dn_spinor_field[0], dn0, VOLUME/2);
    }

    temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);
    Ener[0] = temp;

    temp = square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);
    Ener[0] += temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      ij=0;
      printf("PHMC: Here comes the computation of H_new with \n \n");
      printf("PHMC: At j=%d  P+HMC Final Energy %e \n", ij, mnl->energy1+Ener[0]);
      printf("PHMC: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    mnl->energy1 += Ener[0];
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[0], mnl->energy1);
    }
  } 
  else if(phmc_exact_poly == 1 && g_epsbar == 0.0) {
    phmc_Cpol = mnl->MDPolyLocNormConst;
    for(j = 1; j < (phmc_dop_n_cheby); j++) {
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_min_cconst_nrm(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    mnl->MDPolyRoots[j-1]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);

    temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);
    Ener[0] = temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: Here comes the computation of H_new with \n \n");
      printf("PHMC: At j=%d  P+HMC Final Energy %e \n", ij, mnl->energy1+Ener[0]);
      printf("PHMC: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    mnl->energy1 += Ener[0];
    if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
      printf("PHMC: At j = %d  P=%e +HMC Final Energy %e \n\n", ij, Ener[0], mnl->energy1);
    }
  }

  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called ndpoly_acc for id %d %d dH = %1.4e\n", id, g_running_phmc, mnl->energy1 - mnl->energy0);
  }
  /* END IF PHMC */
  return(mnl->energy1 - mnl->energy0);
}


int init_nd_poly_monomial(const int id) {
  monomial * mnl = &monomial_list[id];
  int j, k, errcode;
  FILE * ifs;
  char title[100];

  phmc_invmaxev = 1.0;
  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  g_kappa = mnl->kappa;
  boundary(g_kappa);
  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
  }
  else {
    phmc_Cpol = mnl->MDPolyLocNormConst;
  }


  /* This is the epsilon parameter */
  mnl->EVMin = mnl->StildeMin / mnl->StildeMax;
  
  /* In the following there is the  "sqrt"  since the value refers to 
     the hermitian Dirac operator (used in EV-computation), namely 
     S = Q Q^dag         
     When  "S"  is applied, we call  phmc_invmaxev  twice !!! */
  if(g_epsbar!=0.0 || phmc_exact_poly==0) mnl->EVMaxInv = 1./(sqrt(mnl->StildeMax));
  else if(g_epsbar==0.0 && phmc_exact_poly==1) mnl->EVMaxInv = 1./mnl->StildeMax;
  phmc_cheb_evmin = mnl->EVMin;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmax = 1.0;

  /* Here we prepare the less precise polynomial first   */
  /* the routine determines a value for phmc_dop_n_cheby */
  degree_of_polynomial_nd(mnl->MDPolyDegree);
  if((g_proc_id == 0) && (g_debug_level > 1)) {
    printf("# monomial %s approximation interval [stilde_min, stilde_max] = [%e, %e]\n", 
	   mnl->name, mnl->StildeMin, mnl->StildeMax);
    printf("# monomial %s degree for P = %d, epsilont = %e, normalisation = %e", 
	   mnl->name, phmc_dop_n_cheby-1, mnl->EVMin, mnl->EVMaxInv);
  }

  /* Chi`s-spinors  memory allocation */
  j = init_chi_spinor_field(VOLUMEPLUSRAND/2, (phmc_dop_n_cheby+1));
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for PHMC Chi fields! Aborting...\n");
    exit(0);
  }

  /* End memory allocation */
  /* Here we prepare the precise polynomial */
  degree_of_Ptilde();

  /* THIS IS THE OVERALL CONSTANT */
  /* write phmc_Cpol as the result of the simple-program files (BigC^(1/2))^1/2 
     since  BigC^(1/2)  is the constant appearing in each factor of the 
     multiplication defining the monomial basis representation of the 
     polinomial in s,  while its square phmc_root  (BigC^(1/2))^1/2  is the 
     constant appearing in the multiplication representing the 
     polinomial in  sqrt(s) .
  */
  if(mnl->MDPolyLocNormConst < 0.0){
    fprintf(stderr, "Error, please specify MDPolyLocNormConst in the input file! Aborting...\n");
#ifdef MPI
    MPI_Finalize();
#endif
    exit(6);
  } 

  mnl->MDPolyRoots = calloc((2*phmc_dop_n_cheby-2),sizeof(_Complex double));

  if((ifs = fopen(mnl->MDPolyRootsFile, "r")) != (FILE*)NULL) {
    if (fgets(title, 100, ifs) == NULL) {
      fprintf(stderr, "Error in reading %s! Aborting...\n", mnl->MDPolyRootsFile);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(6);
    }
    
    /* Here we read in the 2n roots needed for the polinomial in sqrt(s) */
    double *phmc_darray = (double*)mnl->MDPolyRoots;
    for(j = 0; j< 2 * phmc_dop_n_cheby - 2; ++j) {
      errcode = fscanf(ifs, " %d %lf %lf \n", &k, &phmc_darray[2 * j], &phmc_darray[2 * j + 1]);
    }
    fclose(ifs);
  }
  else {
    fprintf(stderr, "File %s is missing! Aborting...\n", mnl->MDPolyRootsFile);
#ifdef MPI
    MPI_Finalize();
#endif
    exit(6);
  }
  
  return(0);
}

void ndpoly_set_global_parameter(monomial * const mnl) {

  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  g_kappa = mnl->kappa;
  boundary(g_kappa);

  phmc_root = mnl->MDPolyRoots;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmin = mnl->EVMin;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmax = 1.0;
 
  return;
}
