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
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/solver.h"
#include "deriv_Sb.h"
#include "operator/tm_operators.h"
#include "chebyshev_polynomial.h"
#include "operator/tm_operators_nd.h"
#include "operator/Hopping_Matrix.h"
#include "phmc.h"
#include "operator/tm_operators_nd.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "reweighting_factor_nd.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "phmc.h"
#include "init/init_chi_spinor_field.h"
#include "solver/matrix_mult_typedef_nd.h"
#include "operator/clover_leaf.h"
#include "operator/clovertm_operators.h"
#include "ndpoly_monomial.h"

extern int phmc_exact_poly;

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void ndpoly_derivative(const int id, hamiltonian_field_t * const hf) {
  double atime, etime;
  int j, k;
  monomial * mnl = &monomial_list[id];
  atime = gettime();
  /* This factor 2 a missing factor 2 in trace_lambda */
  ndpoly_set_global_parameter(mnl, phmc_exact_poly);
  mnl->forcefactor = -phmc_Cpol*mnl->EVMaxInv;
  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */

  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    /* in  g_chi_up_spinor_field[0] (g_chi_dn_spinor_field[0] we expect */
    /* to find the phi field, the pseudo fermion field                  */
    /* i.e. must be equal to mnl->pf (mnl->pf2)                         */

    assign(g_chi_up_spinor_field[0], mnl->pf, VOLUME/2);
    assign(g_chi_dn_spinor_field[0], mnl->pf2, VOLUME/2);

    for(k = 1; k < (mnl->MDPolyDegree-1); k++) {
      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[k], g_chi_dn_spinor_field[k], 
			     g_chi_up_spinor_field[k-1], g_chi_dn_spinor_field[k-1], 
			     mnl->MDPolyRoots[k-1], phmc_Cpol, phmc_invmaxev);
    }
    
    /* Here comes the remaining fields  chi_k ; k=n,...,2n-1  */
    /*They are evaluated step-by-step overwriting the same field (mnl->MDPolyDegree)*/
    
    assign(g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_up_spinor_field[mnl->MDPolyDegree-2], VOLUME/2);
    assign(g_chi_dn_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree-2], VOLUME/2);
    
    for(j=(mnl->MDPolyDegree-1); j>=1; j--) {
      assign(g_chi_up_spinor_field[mnl->MDPolyDegree-1], g_chi_up_spinor_field[mnl->MDPolyDegree], VOLUME/2);
      assign(g_chi_dn_spinor_field[mnl->MDPolyDegree-1], g_chi_dn_spinor_field[mnl->MDPolyDegree], VOLUME/2);
      
      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree], 
			     g_chi_up_spinor_field[mnl->MDPolyDegree-1], g_chi_dn_spinor_field[mnl->MDPolyDegree-1], 
			     mnl->MDPolyRoots[2*mnl->MDPolyDegree-3-j], phmc_Cpol, phmc_invmaxev);
      
      /* Get the even parts of the  (j-1)th  chi_spinors */
      H_eo_tm_ndpsi(mnl->w_fields[0], mnl->w_fields[1], 
		    g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], EO);
      
      /* \delta M_eo sandwitched by  chi[j-1]_e^\dagger  and  chi[2N-j]_o */
      deriv_Sb(EO, mnl->w_fields[0], g_chi_up_spinor_field[mnl->MDPolyDegree], hf, mnl->forcefactor);/* UP */
      deriv_Sb(EO, mnl->w_fields[1], g_chi_dn_spinor_field[mnl->MDPolyDegree], hf, mnl->forcefactor);/* DN */
      
      /* Get the even parts of the  (2N-j)-th  chi_spinors */
      H_eo_tm_ndpsi(mnl->w_fields[0], mnl->w_fields[1], 
		    g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree], EO);
      
      /* \delta M_oe sandwitched by  chi[j-1]_o^\dagger  and  chi[2N-j]_e */
      deriv_Sb(OE, g_chi_up_spinor_field[j-1], mnl->w_fields[0], hf, mnl->forcefactor);
      deriv_Sb(OE, g_chi_dn_spinor_field[j-1], mnl->w_fields[1], hf, mnl->forcefactor);
    }
  } 
  else if(g_epsbar == 0.0) {
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    assign(g_chi_up_spinor_field[0], mnl->pf, VOLUME/2);
    for(k = 1; k < (mnl->MDPolyDegree-1); k++) {
      Qtm_pm_sub_const_nrm_psi(g_chi_up_spinor_field[k],
			       g_chi_up_spinor_field[k-1], 
			       mnl->MDPolyRoots[k-1]);
    }
    assign(g_chi_up_spinor_field[mnl->MDPolyDegree],
	   g_chi_up_spinor_field[mnl->MDPolyDegree-2], VOLUME/2);
    
    for(j = (mnl->MDPolyDegree-1); j >= 1; j--) {
      assign(g_chi_up_spinor_field[mnl->MDPolyDegree-1],
	     g_chi_up_spinor_field[mnl->MDPolyDegree], VOLUME/2);
      
      Qtm_pm_sub_const_nrm_psi(g_chi_up_spinor_field[mnl->MDPolyDegree], 
			       g_chi_up_spinor_field[mnl->MDPolyDegree-1],
			       mnl->MDPolyRoots[2*mnl->MDPolyDegree-3-j]);
      
      Qtm_minus_psi(mnl->w_fields[3],g_chi_up_spinor_field[j-1]); 
      
      H_eo_tm_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[phmc_dop_n_cheby], EO, -1.);
      deriv_Sb(OE, mnl->w_fields[3], mnl->w_fields[2], hf, mnl->forcefactor); 
      
      H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[3], EO, 1.); 
      deriv_Sb(EO, mnl->w_fields[2], g_chi_up_spinor_field[phmc_dop_n_cheby], hf, mnl->forcefactor);
      
      Qtm_minus_psi(mnl->w_fields[3],g_chi_up_spinor_field[mnl->MDPolyDegree]); 

      H_eo_tm_inv_psi(mnl->w_fields[2],mnl->w_fields[3], EO, +1.);
      deriv_Sb(OE, g_chi_up_spinor_field[j-1] , mnl->w_fields[2], hf, mnl->forcefactor); 
      
      H_eo_tm_inv_psi(mnl->w_fields[2], g_chi_up_spinor_field[j-1], EO, -1.); 
      deriv_Sb(EO, mnl->w_fields[2], mnl->w_fields[3], hf, mnl->forcefactor);
    }
  }
  /*
    Normalisation by the largest  EW  is done in update_momenta
    using mnl->forcefactor
  */ 
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}


void ndpoly_heatbath(const int id, hamiltonian_field_t * const hf) {
  int j;
  monomial * mnl = &monomial_list[id];

  ndpoly_set_global_parameter(mnl, phmc_exact_poly);

  // we measure before trajectory!
  if((mnl->rec_ev != 0) && (hf->traj_counter%mnl->rec_ev == 0)) {
    phmc_compute_ev(hf->traj_counter-1, id, &Qtm_pm_ndbipsi);
  }

  mnl->energy0 = 0.;
  random_spinor_field_eo(g_chi_up_spinor_field[0], mnl->rngrepro, RN_GAUSS);
  mnl->energy0 = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);

  if(g_epsbar!=0.0 || phmc_exact_poly == 0) {
    random_spinor_field_eo(g_chi_dn_spinor_field[0], mnl->rngrepro, RN_GAUSS);
    mnl->energy0 += square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);
  } 
  else {
    zero_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2);
  }

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 5)) {
    printf("# NDPOLY: OLD Energy  DN + UP %e \n\n", mnl->energy0);
  }

  if(phmc_exact_poly==0){
    Qtm_ndpsi(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
		    g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0]);
 
    for(j = 1; j < (mnl->MDPolyDegree); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], 
			mnl->MDPolyRoots[mnl->MDPolyDegree-2+j], phmc_Cpol, phmc_invmaxev);
    }
    Ptilde_ndpsi(g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], mnl->PtildeCoefs, 
		 mnl->PtildeDegree, g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], &Qtm_pm_ndpsi);
  } 
  else if( phmc_exact_poly==1 && g_epsbar!=0.0) {
    /* Attention this is Q * tau1, up/dn are exchanged in the input spinor  */
    /* this is used as an preconditioner */
    Qtm_ndpsi(g_chi_up_spinor_field[1],g_chi_dn_spinor_field[1],
		    g_chi_dn_spinor_field[0],g_chi_up_spinor_field[0]);

    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);

    /* solve Q*tau1*P(Q^2) *x=y */
    cg_her_nd(g_chi_up_spinor_field[1],g_chi_dn_spinor_field[1],
	      g_chi_up_spinor_field[0],g_chi_dn_spinor_field[0],
	      1000,1.e-16,0,VOLUME/2, Qtau1_P_ndpsi);

    /*  phi= Bdagger phi  */
    for(j = 1; j < (mnl->MDPolyDegree); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
      Q_tau1_sub_const_ndpsi(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1],
			g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0],
			mnl->MDPolyRoots[mnl->MDPolyDegree-2+j], phmc_Cpol, phmc_invmaxev);
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
    for(j = 1; j < (mnl->MDPolyDegree); j++){
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_sub_const_nrm_psi(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    mnl->MDPolyRoots[mnl->MDPolyDegree-2+j]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
  }

  assign(mnl->pf, g_chi_up_spinor_field[0], VOLUME/2);
  assign(mnl->pf2, g_chi_dn_spinor_field[0], VOLUME/2);

  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called ndpoly_heatbath for id %d \n", id);
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

  ndpoly_set_global_parameter(mnl, phmc_exact_poly);
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
    for(j = 1; j <= (mnl->MDPolyDegree-1); j++) {
      /* Change this name !!*/
      Q_tau1_sub_const_ndpsi(up1, dn1, up0, dn0, mnl->MDPolyRoots[j-1], phmc_Cpol, phmc_invmaxev);

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

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 4)) {
      printf("# NDPOLY: At j=%d H before H-correction %e \n", ij, Ener[ij]);
    }
    
    /* Here comes the loop for the evaluation of A, A^2, ...  */
    for(j = 1; j < 8; j++){ /* To omit corrections just set  j<1 */
      
      if(j % 2){ /*  Chi[j] = ( Qdag P  Ptilde ) Chi[j-1]  */ 
	Ptilde_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		      mnl->PtildeCoefs, mnl->PtildeDegree, 
		     g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], &Qtm_pm_ndpsi);
	Ptilde_ndpsi(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		     mnl->MDPolyCoefs, mnl->MDPolyDegree, 
		     g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], &Qtm_pm_ndpsi);

	Qtm_dagger_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
      }
      else { /*  Chi[j] = ( Ptilde P Q ) Chi[j-1]  */ 
	Qtm_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
			g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
	Ptilde_ndpsi(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		     mnl->MDPolyCoefs, mnl->MDPolyDegree, g_chi_up_spinor_field[j], 
		     g_chi_dn_spinor_field[j], &Qtm_pm_ndpsi);
	Ptilde_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		      mnl->PtildeCoefs, mnl->PtildeDegree, 
		      g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], &Qtm_pm_ndpsi);
      }

      Ener[j] = Ener[j-1] + Ener[0];
      sgn = -1.0;
      for(ij = 1; ij < j; ij++){
	fact = factor[j] / (factor[ij] * factor[j-ij]);
	if((g_proc_id == g_stdio_proc) && (g_debug_level > 4)) {
	  printf("# NDPOLY: Here  j=%d  and  ij=%d   sign=%f  fact=%f \n", j ,ij, sgn, fact);
	}
	Ener[j] += sgn*fact*Ener[ij];
	sgn = -sgn;
      }
      temp = square_norm(g_chi_up_spinor_field[j], VOLUME/2, 1);
      temp += square_norm(g_chi_dn_spinor_field[j], VOLUME/2, 1);
      if((g_proc_id == g_stdio_proc) && (g_debug_level > 4)) {
	printf("# NDPOLY: Here  j=%d   sign=%f  temp=%e \n", j, sgn, temp);
      }

      Ener[j] += sgn*temp;

      Diff = fabs(Ener[j] - Ener[j-1]);
      if((g_proc_id == g_stdio_proc) && (g_debug_level > 0)) {
	printf("# NDPOLY: H-Correction after %d steps: %e \n", j, Diff);
      }

      if(Diff < mnl->PrecisionHfinal) {
	break;
      }
    }
    mnl->energy1 += Ener[ij];  /* this is quite sticky */
  } 
  else if(phmc_exact_poly==1 && g_epsbar!=0.0) {
    /* B(Q*tau1) */
    for(j = 1; j <= (mnl->MDPolyDegree-1); j++){
      Q_tau1_sub_const_ndpsi(up1, dn1, up0, dn0, mnl->MDPolyRoots[j-1], phmc_Cpol, phmc_invmaxev);

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

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 4)) {
      ij=0;
      printf("# NDPOLY: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    mnl->energy1 += Ener[0];
  } 
  else if(phmc_exact_poly == 1 && g_epsbar == 0.0) {
    for(j = 1; j < (mnl->MDPolyDegree); j++) {
      assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
      Qtm_pm_sub_const_nrm_psi(g_chi_up_spinor_field[1],
			    g_chi_up_spinor_field[0],
			    mnl->MDPolyRoots[j-1]);
    }
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);

    temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);
    Ener[0] = temp;

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 4)) {
      printf("# NDPOLY: At j=%d  PHMC Only Final Energy %e \n", ij, Ener[0]);
    }

    mnl->energy1 += Ener[0];
  }

  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called ndpoly_acc for id %d %d dH = %1.10e\n", id, g_running_phmc, mnl->energy1 - mnl->energy0);
  }
  /* END IF PHMC */
  return(mnl->energy1 - mnl->energy0);
}


int init_ndpoly_monomial(const int id) {
  monomial * mnl = &monomial_list[id];
  int j, k, errcode;
  FILE * ifs;
  double *phmc_darray;
  char title[100];
  matrix_mult_nd Qsq = &Qtm_pm_ndpsi;
  double atime, etime;

  atime = gettime();
  if(mnl->type == NDCLOVER) {
    Qsq = &Qsw_pm_ndpsi;
    init_sw_fields();
    sw_term((const su3 **)g_gauge_field, mnl->kappa, mnl->c_sw); 
    sw_invert_nd(mnl->mubar*mnl->mubar - mnl->epsbar*mnl->epsbar);
  }

  phmc_invmaxev = 1.0;
  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  g_kappa = mnl->kappa;
  g_c_sw = mnl->c_sw;
  boundary(g_kappa);
  if (g_epsbar!=0.0 || phmc_exact_poly==0){
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
  }
  else {
    phmc_Cpol = mnl->MDPolyLocNormConst;
  }

  /* This is the epsilon parameter */
  mnl->EVMin = mnl->StildeMin / mnl->StildeMax;
  mnl->EVMax = 1.;
  /* In the following there is the  "sqrt"  since the value refers to 
     the hermitian Dirac operator (used in EV-computation), namely 
     S = Q Q^dag         
     When  "S"  is applied, we call  phmc_invmaxev  twice !!! */
  if(g_epsbar!=0.0 || phmc_exact_poly==0) mnl->EVMaxInv = 1./(sqrt(mnl->StildeMax));
  else if(g_epsbar==0.0 && phmc_exact_poly==1) mnl->EVMaxInv = 1./mnl->StildeMax;
  phmc_cheb_evmin = mnl->EVMin;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmax = 1.0;

  /* Here we prepare the less precise MD polynomial first   */
  degree_of_polynomial_nd(&mnl->MDPolyDegree, &mnl->MDPolyCoefs,
			  mnl->EVMin, mnl->EVMax,
			  Qsq, mnl->rngrepro);
  phmc_dop_n_cheby = mnl->MDPolyDegree;
  phmc_dop_cheby_coef = mnl->MDPolyCoefs;
  if((g_proc_id == 0) && (g_debug_level > 1)) {
    printf("# monomial %s approximation interval [stilde_min, stilde_max] = [%e, %e]\n", 
	   mnl->name, mnl->StildeMin, mnl->StildeMax);
    printf("# monomial %s degree for P = %d, epsilont = %e, normalisation = %e", 
	   mnl->name, mnl->MDPolyDegree-1, mnl->EVMin, mnl->EVMaxInv);
  }

  /* Chi`s-spinors  memory allocation */
  j = init_chi_spinor_field(VOLUMEPLUSRAND/2, (mnl->MDPolyDegree+1));
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for PHMC Chi fields! Aborting...\n");
    exit(0);
  }

  /* End memory allocation */
  /* Here we prepare the precise polynomial Ptilde */
  degree_of_Ptilde(&mnl->PtildeDegree, &mnl->PtildeCoefs, 
		   mnl->EVMin, mnl->EVMax, mnl->MDPolyDegree, 
		   mnl->PrecisionPtilde, Qsq, mnl->rngrepro);
  phmc_ptilde_cheby_coef = mnl->PtildeCoefs;
  phmc_ptilde_n_cheby = mnl->PtildeDegree;

  /* THIS IS THE OVERALL CONSTANT */
  /* write phmc_Cpol as the result of the simple-program files (BigC^(1/2))^1/2 
     since  BigC^(1/2)  is the constant appearing in each factor of the 
     multiplication defining the monomial basis representation of the 
     polinomial in s,  while its square phmc_root  (BigC^(1/2))^1/2  is the 
     constant appearing in the multiplication representing the 
     polinomial in  sqrt(s) .
  */
  if(mnl->MDPolyLocNormConst < 0.0){
    fprintf(stderr, "Error, please specify LocNormConst in the input file! Aborting...\n");
#ifdef MPI
    MPI_Finalize();
#endif
    exit(6);
  } 

  mnl->MDPolyRoots = calloc((2*mnl->MDPolyDegree-2),sizeof(_Complex double));

  if((ifs = fopen(mnl->MDPolyRootsFile, "r")) != (FILE*)NULL) {
    if (fgets(title, 100, ifs) == NULL) {
      fprintf(stderr, "Error in reading %s! Aborting...\n", mnl->MDPolyRootsFile);
#ifdef MPI
      MPI_Finalize();
#endif
      exit(6);
    }
    
    /* Here we read in the 2n roots needed for the polinomial in sqrt(s) */
    phmc_darray = (double*)mnl->MDPolyRoots;
    for(j = 0; j< 2 * mnl->MDPolyDegree - 2; ++j) {
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
  etime = gettime();
  if(g_debug_level > 0 && g_proc_id == 0) {
    printf("# Time for init %s monomial: %e s\n", mnl->name, etime-atime);
  }
  return(0);
}

void ndpoly_set_global_parameter(monomial * const mnl, const int exact) {

  g_mubar = mnl->mubar;
  g_epsbar = mnl->epsbar;
  g_kappa = mnl->kappa;
  g_c_sw = mnl->c_sw;
  boundary(g_kappa);

  if (g_epsbar!=0.0 || exact == 0){
    phmc_Cpol = sqrt(mnl->MDPolyLocNormConst);
  }
  else {
    phmc_Cpol = mnl->MDPolyLocNormConst;
  }

  phmc_root = mnl->MDPolyRoots;
  phmc_cheb_evmin = mnl->EVMin;
  phmc_invmaxev = mnl->EVMaxInv;
  phmc_cheb_evmax = 1.0;

  phmc_dop_n_cheby = mnl->MDPolyDegree;
  phmc_dop_cheby_coef = mnl->MDPolyCoefs;

  phmc_ptilde_cheby_coef = mnl->PtildeCoefs;
  phmc_ptilde_n_cheby = mnl->PtildeDegree;
 
  return;
}
