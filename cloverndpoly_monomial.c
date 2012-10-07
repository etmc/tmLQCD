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
#include "tm_operators_nd.h"
#include "Hopping_Matrix.h"
#include "phmc.h"
#include "tm_operators_nd.h"
#include "chebyshev_polynomial_nd.h"
#include "Ptilde_nd.h"
#include "reweighting_factor_nd.h"
#include "monomial.h"
#include "hamiltonian_field.h"
#include "boundary.h"
#include "phmc.h"
#include "init_chi_spinor_field.h"
#include "clover_leaf.h"
#include "cloverndpoly_monomial.h"

/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

void cloverndpoly_derivative(const int id, hamiltonian_field_t * const hf) {
  int j, k;
  monomial * mnl = &monomial_list[id];

  for(int i = 0; i < VOLUME; i++) { 
    for(int mu = 0; mu < 4; mu++) { 
      _su3_zero(swm[i][mu]);
      _su3_zero(swp[i][mu]);
    }
  }
  ndpoly_set_global_parameter(mnl, 0);
  
  // we compute the clover term (1 + T_ee(oo)) for all sites x
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  // we invert it for the even sites only
  sw_invert_nd(EE, mnl->mu);


  /* This factor 2 a missing factor 2 in trace_lambda */
  ndpoly_set_global_parameter(mnl, 0);
  mnl->forcefactor = -phmc_Cpol*mnl->EVMaxInv;

  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */

  /* Here comes the definitions for the chi_j fields */
  /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
  /* in  g_chi_up_spinor_field[0] (g_chi_dn_spinor_field[0] we expect */
  /* to find the phi field, the pseudo fermion field                  */
  /* i.e. must be equal to mnl->pf (mnl->pf2)                         */
  
  assign(g_chi_up_spinor_field[0], mnl->pf, VOLUME/2);
  assign(g_chi_dn_spinor_field[0], mnl->pf2, VOLUME/2);
  
  for(k = 1; k < (mnl->MDPolyDegree-1); k++) {
    Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[k], g_chi_dn_spinor_field[k], 
			     g_chi_up_spinor_field[k-1], g_chi_dn_spinor_field[k-1], 
			     mnl->MDPolyRoots[k-1]);
  }
  
  /* Here comes the remaining fields  chi_k ; k=n,...,2n-1  */
  /*They are evaluated step-by-step overwriting the same field (mnl->MDPolyDegree)*/
  
  assign(g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_up_spinor_field[mnl->MDPolyDegree-2], VOLUME/2);
  assign(g_chi_dn_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree-2], VOLUME/2);
  
  for(j=(mnl->MDPolyDegree-1); j>=1; j--) {
    assign(g_chi_up_spinor_field[mnl->MDPolyDegree-1], g_chi_up_spinor_field[mnl->MDPolyDegree], VOLUME/2);
    assign(g_chi_dn_spinor_field[mnl->MDPolyDegree-1], g_chi_dn_spinor_field[mnl->MDPolyDegree], VOLUME/2);
    
    Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree], 
			     g_chi_up_spinor_field[mnl->MDPolyDegree-1], g_chi_dn_spinor_field[mnl->MDPolyDegree-1], 
			     mnl->MDPolyRoots[2*mnl->MDPolyDegree-3-j]);
    
    /* Get the even parts of the  (j-1)th  chi_spinors */
    H_eo_sw_ndpsi(mnl->w_fields[0], mnl->w_fields[1], 
		  g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], EO);
    
    /* \delta M_eo sandwitched by  chi[j-1]_e^\dagger  and  chi[2N-j]_o */
    deriv_Sb(EO, mnl->w_fields[0], g_chi_up_spinor_field[mnl->MDPolyDegree], hf, mnl->forcefactor);      /* UP */
    deriv_Sb(EO, mnl->w_fields[1], g_chi_dn_spinor_field[mnl->MDPolyDegree], hf, mnl->forcefactor);    /* DN */
    
    /* Get the even parts of the  (2N-j)-th  chi_spinors */
    H_eo_sw_ndpsi(mnl->w_fields[0], mnl->w_fields[1], 
	    g_chi_up_spinor_field[mnl->MDPolyDegree], g_chi_dn_spinor_field[mnl->MDPolyDegree], EO);
    
    /* \delta M_oe sandwitched by  chi[j-1]_o^\dagger  and  chi[2N-j]_e */
    deriv_Sb(OE, g_chi_up_spinor_field[j-1], mnl->w_fields[0], hf, mnl->forcefactor);
    deriv_Sb(OE, g_chi_dn_spinor_field[j-1], mnl->w_fields[1], hf, mnl->forcefactor);
  }
  return;
}


void cloverndpoly_heatbath(const int id, hamiltonian_field_t * const hf) {
  int j;
  double temp;
  monomial * mnl = &monomial_list[id];

  ndpoly_set_global_parameter(mnl, 0);
  g_mu3 = 0.;
  init_sw_fields();
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert_nd(EE, mnl->mu);

  mnl->energy0 = 0.;
  random_spinor_field(g_chi_up_spinor_field[0], VOLUME/2, mnl->rngrepro);
  mnl->energy0 = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);

  random_spinor_field(g_chi_dn_spinor_field[0], VOLUME/2, mnl->rngrepro);
  mnl->energy0 += square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Here comes the computation of H_old with \n \n");
    printf("PHMC: First: random spinors and their norm  \n ");
    printf("PHMC: OLD Ennergy UP %e \n", mnl->energy0);
    printf("PHMC: OLD Energy  DN + UP %e \n\n", mnl->energy0);
  }

  Qsw_ndpsi(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
		  g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0]);
  
  for(j = 1; j < (mnl->MDPolyDegree); j++){
    assign(g_chi_up_spinor_field[0], g_chi_up_spinor_field[1], VOLUME/2);
    assign(g_chi_dn_spinor_field[0], g_chi_dn_spinor_field[1], VOLUME/2);
    
    Qsw_tau1_sub_const_ndpsi(g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1], 
			 g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], 
			 mnl->MDPolyRoots[mnl->MDPolyDegree-2+j]);
  }
  Ptilde_ndpsi(g_chi_up_spinor_field[0], g_chi_dn_spinor_field[0], mnl->PtildeCoefs, 
	       mnl->PtildeDegree, g_chi_up_spinor_field[1], g_chi_dn_spinor_field[1]);
  
  assign(mnl->pf, g_chi_up_spinor_field[0], VOLUME/2);
  assign(mnl->pf2, g_chi_dn_spinor_field[0], VOLUME/2);
  
  temp = square_norm(g_chi_up_spinor_field[0], VOLUME/2, 1);
  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)) {
    printf("PHMC: Then: evaluate Norm of pseudofermion heatbath BHB \n ");
    printf("PHMC: Norm of BHB up squared %e \n", temp);
  }

  temp += square_norm(g_chi_dn_spinor_field[0], VOLUME/2, 1);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 2)){
    printf("PHMC: Norm of BHB up + BHB dn squared %e \n\n", temp);
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverndpoly_heatbath for id %d with g_running_phmc = %d\n", id, g_running_phmc);
  }
  return;
}


double cloverndpoly_acc(const int id, hamiltonian_field_t * const hf) {
  int j, ij=0;
  double temp, sgn, fact, Diff;
  double Ener[8];
  double factor[8];
  monomial * mnl = &monomial_list[id];
  spinor *up0, *dn0, *up1, *dn1, *dummy;

  ndpoly_set_global_parameter(mnl, 0);
  g_mu3 = 0.;
  init_sw_fields();
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  sw_invert_nd(EE, mnl->mu);

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

  for(j = 1; j <= (mnl->MDPolyDegree-1); j++) {
    /* Change this name !!*/
    Qsw_tau1_sub_const_ndpsi(up1, dn1, up0, dn0, mnl->MDPolyRoots[j-1]);
    
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
      Ptilde_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		   mnl->PtildeCoefs, mnl->PtildeDegree, 
		   g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
      Ptilde_ndpsi(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		   mnl->MDPolyCoefs, mnl->MDPolyDegree, 
		   g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j]);
      Qsw_dagger_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		       g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
    }
    else { /*  Chi[j] = ( Ptilde P Q ) Chi[j-1]  */ 
      Qsw_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1]);
      Ptilde_ndpsi(g_chi_up_spinor_field[j-1], g_chi_dn_spinor_field[j-1], 
		   mnl->MDPolyCoefs, mnl->MDPolyDegree, g_chi_up_spinor_field[j], 
		   g_chi_dn_spinor_field[j]);
      Ptilde_ndpsi(g_chi_up_spinor_field[j], g_chi_dn_spinor_field[j], 
		   mnl->PtildeCoefs, mnl->PtildeDegree, 
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
    
    if(Diff < mnl->PrecisionHfinal) {
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
  
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called cloverndpoly_acc for id %d %d dH = %1.4e\n", id, g_running_phmc, mnl->energy1 - mnl->energy0);
  }
  return(mnl->energy1 - mnl->energy0);
}


