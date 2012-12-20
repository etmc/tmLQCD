/***********************************************************************
 * 
 * Copyright (C) 2010 Andreas Nube
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
#include "global.h"
#include "start.h"
#include "gettime.h"
#include "read_input.h"
#include "monomial/monomial.h"
#include "poly_monomial.h"
#include "boundary.h"
#include "linalg/square_norm.h"
#include "linalg/assign.h"
#include "linalg/mul_r.h"
#include "linalg/diff.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "operator/tm_operators.h"
#include "solver/solver.h"
#include "solver/chrono_guess.h"
#include "solver/eigenvalues.h"
#include "operator/tm_operators_nd.h"
#include "operator/Hopping_Matrix.h"
#include "hamiltonian_field.h"
#include "phmc.h"



inline void setPhmcVars(monomial *mnl){
  phmc_invmaxev=1.0/mnl->MDPolyLmax;
  phmc_dop_n_cheby=(mnl->MDPolyDegree/2)+1;
  phmc_Cpol=mnl->MDPolyLocNormConst;
  phmc_root=mnl->MDPolyRoots;
}

void poly_derivative(const int id, hamiltonian_field_t * const hf){
  double atime, etime;
  monomial * mnl = &monomial_list[id];
  int k,j;
  int degreehalf=mnl->MDPolyDegree/2;

  spinor** chi_spinor_field=mnl->MDPoly_chi_spinor_fields;

  atime = gettime();
  (*mnl).forcefactor = -mnl->MDPolyLocNormConst/mnl->MDPolyLmax;


  /* push and set phmc vars */
  pushPhmcVars();
  setPhmcVars(mnl);


  if(mnl->even_odd_flag){



    if(mnl->MDPolyDetRatio==1){
      g_mu=mnl->mu2;
      boundary(mnl->kappa2);
      Qtm_plus_psi(chi_spinor_field[0],mnl->pf);
    } else {
      assign(chi_spinor_field[0],mnl->pf,VOLUME/2);
    }


    g_mu=mnl->mu;
    boundary(mnl->kappa);

    
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    for(k = 0; k < degreehalf-1 ; k++) {
      Qtm_pm_sub_const_nrm_psi(chi_spinor_field[k+1],
				 chi_spinor_field[k], 
				 mnl->MDPolyRoots[k]);
    }



    assign(chi_spinor_field[degreehalf+1],
	   chi_spinor_field[degreehalf-1], VOLUME/2);
  
    /* loop over monoms */
    for(j=degreehalf; j>=1; j--) {

      assign(chi_spinor_field[degreehalf],
	     chi_spinor_field[degreehalf+1], VOLUME/2);
      
      Qtm_pm_sub_const_nrm_psi(chi_spinor_field[degreehalf+1], 
			    chi_spinor_field[degreehalf],
			    mnl->MDPolyRoots[mnl->MDPolyDegree-(j+1)]);
      

      Qtm_minus_psi(mnl->w_fields[1],chi_spinor_field[j-1]); 

      H_eo_tm_inv_psi(mnl->w_fields[0], chi_spinor_field[degreehalf+1], EO, -1.);
      deriv_Sb(OE, mnl->w_fields[1], mnl->w_fields[0], hf, mnl->forcefactor); 
      
      H_eo_tm_inv_psi(mnl->w_fields[0], mnl->w_fields[1], EO, 1.); 
      deriv_Sb(EO, mnl->w_fields[0], chi_spinor_field[degreehalf+1], hf, mnl->forcefactor);
    
      Qtm_minus_psi(mnl->w_fields[1],chi_spinor_field[degreehalf+1]); 

      H_eo_tm_inv_psi(mnl->w_fields[0],mnl->w_fields[1], EO, +1.);
      deriv_Sb(OE, chi_spinor_field[j-1] , mnl->w_fields[0], hf, mnl->forcefactor); 
      
      H_eo_tm_inv_psi(mnl->w_fields[0], chi_spinor_field[j-1], EO, -1.); 
      deriv_Sb(EO, mnl->w_fields[0], mnl->w_fields[1], hf, mnl->forcefactor);
    }


    if(mnl->MDPolyDetRatio==1){
      /******************************************
      * multiply with the last missing monomial *
      * such that we get an evaluation of P     *
      ******************************************/
      Qtm_pm_sub_const_nrm_psi(chi_spinor_field[degreehalf], 
			    chi_spinor_field[degreehalf+1],
			    mnl->MDPolyRoots[mnl->MDPolyDegree-1]);
      
    /* devide by this factor cause its multiplied again in update_fermion_momenta see comment below */
    mul_r(chi_spinor_field[degreehalf],
	  1./mnl->MDPolyLocNormConst*mnl->MDPolyLmax,
	  chi_spinor_field[degreehalf],
	  VOLUME/2);


      g_mu=mnl->mu2;
      boundary(mnl->kappa2);

      H_eo_tm_inv_psi(mnl->w_fields[0],chi_spinor_field[degreehalf], EO, -1.);
      deriv_Sb(OE, mnl->pf , mnl->w_fields[0], hf, mnl->forcefactor);

      H_eo_tm_inv_psi(mnl->w_fields[0], mnl->pf, EO, +1.);
      deriv_Sb(EO, mnl->w_fields[0], chi_spinor_field[degreehalf], hf, mnl->forcefactor);
    }
  } 
  else {
    if(g_proc_id == 0) {
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }
    
    g_mu = g_mu1;
    boundary(g_kappa);
    popPhmcVars();
    
    return;
  }

  /* restore all changed global vars */
  g_mu = g_mu1;
  boundary(g_kappa);
  popPhmcVars();
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}

double poly_acc(const int id, hamiltonian_field_t * const hf){

  monomial * mnl = &monomial_list[id];
  int j;
  double diff;
  int no_eigenvalues=-1;
  double atime, etime;
  atime = gettime();
  if(mnl->even_odd_flag) {
    if(mnl->MDPolyDetRatio==1) {
      g_mu = mnl->mu2;
      boundary(mnl->kappa2);
      Qtm_plus_psi(mnl->w_fields[1],mnl->pf);
    } 
    else {
      assign(mnl->w_fields[1],mnl->pf,VOLUME/2);
    }

    g_mu = mnl->mu;
    boundary(mnl->kappa);

    /* push and set phmc vars */
    pushPhmcVars();
    setPhmcVars(mnl);

    /* apply B */
    for(j = 0; j < mnl->MDPolyDegree/2; j++){
      assign(mnl->w_fields[0], mnl->w_fields[1], VOLUME/2);
      Qtm_pm_sub_const_nrm_psi(mnl->w_fields[1],
			    mnl->w_fields[0],
			    mnl->MDPolyRoots[j]);
    }

    mnl->energy1 =  square_norm(mnl->w_fields[1], VOLUME/2,1);

    /* calculate evs */
    if (compute_evs != 0) {
      no_eigenvalues=10;
      eigenvalues(&no_eigenvalues, mnl->maxiter, eigenvalue_precision,
                  0/* compute minimal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);

      no_eigenvalues=1;
      eigenvalues(&no_eigenvalues, mnl->maxiter, eigenvalue_precision,
                  1/* compute maximal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);
    }


    /* restore global phmc vars */
    popPhmcVars();

    
    /* return the energy differnce */
    g_mu = g_mu1;
    boundary(g_kappa);


  
    if(g_proc_id == 0 && g_debug_level > 3) {
      fprintf(stderr," Poly energy1     = %e \n" , mnl->energy1);
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
      diff = mnl->energy1 - mnl->energy0;
      fprintf(stderr," Poly energy diff = %e \n" , diff);
    }
  } 
  else {
    if(g_proc_id == 0) {
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }
    
    g_mu = g_mu1;
    boundary(g_kappa);

    return NAN;
  }
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called poly_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy1 - mnl->energy0);
    }
  }
  return (mnl->energy1 - mnl->energy0);
}

void poly_heatbath(const int id, hamiltonian_field_t * const hf){
  monomial * mnl = &monomial_list[id];
  int j;
  double atime, etime;
  atime = gettime();
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;

  g_mu = mnl->mu;
  boundary(mnl->kappa);
  
  /* push and set phmc vars */
  pushPhmcVars();
  setPhmcVars(mnl);

  if(mnl->even_odd_flag) {


    random_spinor_field_eo(mnl->w_fields[0], mnl->rngrepro, RN_GAUSS);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME/2, 1);

    if(g_proc_id == 0 && g_debug_level > 3) {
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
    }

    /* calculate the phmc hamiltonian */
    Qtm_pm_psi(mnl->w_fields[1], mnl->w_fields[0]);

    /* solve (Q+)*(Q-)*P((Q+)*(Q-)) *x=y */
    cg_her(mnl->w_fields[0], mnl->w_fields[1],
	   1000,mnl->accprec,g_relative_precision_flag,VOLUME/2, Qtm_pm_Ptm_pm_psi);
    
    /*  phi= Bdagger phi  */
    for(j = 0; j < (mnl->MDPolyDegree/2); j++){
      assign(mnl->w_fields[1], mnl->w_fields[0], VOLUME/2);
      Qtm_pm_sub_const_nrm_psi(mnl->w_fields[0],
				 mnl->w_fields[1],
				 mnl->MDPolyRoots[mnl->MDPolyDegree/2+j]);
    }


    if(mnl->MDPolyDetRatio==1){
      g_mu = mnl->mu2;
      boundary(mnl->kappa2);
      zero_spinor_field(mnl->pf,VOLUME/2);
      mnl->iter0 = cg_her(mnl->w_fields[1], mnl->w_fields[0], mnl->maxiter, mnl->accprec, g_relative_precision_flag,
			  VOLUME/2, &Qtm_pm_psi);
      Qtm_minus_psi(mnl->pf, mnl->w_fields[1]);

      chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
			  mnl->csg_N, &mnl->csg_n, VOLUME/2);
      
    } else {
      /* store constructed phi field */
      assign(mnl->pf, mnl->w_fields[0], VOLUME/2);
    }
    
  }
  else {
    /* not implemented */
    fprintf(stderr,"Error: non even odd preconditioned \"light\" phmc not implemented \n");
  }

  g_mu = g_mu1;
  boundary(g_kappa);
  popPhmcVars();

  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called poly_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}
