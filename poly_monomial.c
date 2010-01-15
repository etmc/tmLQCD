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
/*
.
.
.
*/

#include "global.h"
#include "start.h"
#include "read_input.h"
#include "monomial.h"
#include "poly_monomial.h"
#include "boundary.h"
#include "linalg/square_norm.h"
#include "linalg/assign.h"
#include "linalg/mul_r.h"
#include "linalg/diff.h"
#include "linalg_eo.h"
#include "deriv_Sb.h"
#include "linsolve.h"
#include "tm_operators.h"
#include "solver/solver.h"
#include "solver/chrono_guess.h"
#include "solver/eigenvalues.h"
#include "Nondegenerate_Matrix.h"
#include "Hopping_Matrix.h"




monomial *actual_mnl=NULL;
int poly_monomial_V=0;

/* this is neccessary for the calculation of the polynomial */

void Poly_Qtm_pm_min_cconst_nrm(spinor * const l, spinor * const k,
			   const complex z){
  static su3_vector phi1;
  spinor *r,*s;
  int ix;
  Qtm_pm_psi(l,k);
  mul_r(l, 1.0/ actual_mnl->MDPolyLmax, l, poly_monomial_V);

  /*  AND FINALLY WE SUBSTRACT THE C-CONSTANT  */


  /************ loop over all lattice sites ************/
  for(ix = 0; ix < (poly_monomial_V); ix++){

    r=l + ix;
    s=k + ix;

    _complex_times_vector(phi1, z, (*s).s0);
    _vector_sub_assign((*r).s0, phi1);
    _complex_times_vector(phi1, z, (*s).s1);
    _vector_sub_assign((*r).s1, phi1);
    _complex_times_vector(phi1, z, (*s).s2);
    _vector_sub_assign((*r).s2, phi1);
    _complex_times_vector(phi1, z, (*s).s3);
    _vector_sub_assign((*r).s3, phi1);
    
  }

  mul_r(l, actual_mnl->MDPolyLocNormConst, l, poly_monomial_V);
  return;
}



void Poly_Ptm_pm_psi(spinor * const l, spinor * const k){
  
  int j;
  spinor *spinDum;
  spinDum=g_spinor_field[DUM_MATRIX+2];
  
  assign(spinDum,k,poly_monomial_V);
  
  
  for(j=0; j<actual_mnl->MDPolyDegree; j++){
    if(j>0) {
      assign(spinDum,l,poly_monomial_V);
    }
    
    Poly_Qtm_pm_min_cconst_nrm(l,spinDum,actual_mnl->MDPolyRoots[j]);
  }
  return;
}

/* **********************************************
 * Qpm * P(Qpm)
 * this operator is neccessary for the inverter 
 ************************************************/

void Poly_Qtm_pm_Ptm_pm_psi(spinor * const l, spinor * const k){
  spinor * spinDum;
  
  spinDum=g_spinor_field[DUM_MATRIX+3];
  Poly_Ptm_pm_psi(l,k);
  assign(spinDum,l,poly_monomial_V);
  Qtm_pm_psi(l,spinDum);
  return;
}


inline void poly_monomial_c_msg(const char *msg){
  return;
  static int c;
  c=5;
  while(c--)
    printf("        ***********************************\n");
  c=5;
  while(c--)
    printf("%s\n",msg);
  c=5;
  while(c--)
    printf("        ***********************************\n");
}

void poly_derivative(const int id){

  monomial * mnl = &monomial_list[id];
  int k,j;
  int degreehalf=mnl->MDPolyDegree/2;
  poly_monomial_c_msg("Doing Derivative for degenerated Polynomial ...");

  spinor** chi_spinor_field=mnl->MDPoly_chi_spinor_fields;

  g_mu=mnl->mu;

  boundary(mnl->kappa);


  (*mnl).forcefactor = -2.*mnl->MDPolyLocNormConst/mnl->MDPolyLmax;

  actual_mnl=mnl;




  if(mnl->even_odd_flag){

    poly_monomial_V=VOLUME/2;


    assign(chi_spinor_field[0],mnl->pf,VOLUME/2);

    
    /* Here comes the definitions for the chi_j fields */
    /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
    for(k = 0; k < degreehalf-1 ; k++) {
      Poly_Qtm_pm_min_cconst_nrm(chi_spinor_field[k+1],
				 chi_spinor_field[k], 
				 mnl->MDPolyRoots[k]);
    }



    assign(chi_spinor_field[degreehalf+1],
	   chi_spinor_field[degreehalf-1], VOLUME/2);
  
  for(j=degreehalf; j>=1; j--) {
    
/*     if(j==degreehalf){ */
/*       mul_r(chi_spinor_field[degreehalf], -phmc_Cpol*phmc_invmaxev, */
/* 	    chi_spinor_field[degreehalf+1], VOLUME/2); */
/*     } else { */
    assign(chi_spinor_field[degreehalf],
	   chi_spinor_field[degreehalf+1], VOLUME/2);
/*     } */
    
    Poly_Qtm_pm_min_cconst_nrm(chi_spinor_field[degreehalf+1], 
			  chi_spinor_field[degreehalf],
			       mnl->MDPolyRoots[mnl->MDPolyDegree-(j+1)]);
    

    Qtm_minus_psi(g_spinor_field[DUM_DERI+3],chi_spinor_field[j-1]); 
    
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], chi_spinor_field[degreehalf+1], EO, -1.);
    deriv_Sb(OE, g_spinor_field[DUM_DERI+3], g_spinor_field[DUM_DERI+2]); 
    
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], EO, 1.); 
    deriv_Sb(EO, g_spinor_field[DUM_DERI+2], chi_spinor_field[degreehalf+1]);
    
    
   
    
    Qtm_minus_psi(g_spinor_field[DUM_DERI+3],chi_spinor_field[degreehalf+1]); 
    
    
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2],g_spinor_field[DUM_DERI+3], EO, +1.);
    deriv_Sb(OE, chi_spinor_field[j-1] , g_spinor_field[DUM_DERI+2]); 
    
    H_eo_tm_inv_psi(g_spinor_field[DUM_DERI+2], chi_spinor_field[j-1], EO, -1.); 
    deriv_Sb(EO, g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3]);
    
  }




  } else {
    if(g_proc_id==0){
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }

    g_mu = g_mu1;
    boundary(g_kappa);

    return;
  }


  g_mu = g_mu1;
  boundary(g_kappa);


}

double poly_acc(const int id){

  monomial * mnl = &monomial_list[id];
  int j;
  spinor* spinor1=g_spinor_field[2];
  spinor* spinor2=g_spinor_field[3];
  double diff;
  int no_eigenvalues=-1;

  poly_monomial_c_msg("Doing Accept/Reject preparations for degenerate Polynomial");


  g_mu = mnl->mu;
  boundary(mnl->kappa);

  actual_mnl=mnl;

  if(mnl->even_odd_flag){

    poly_monomial_V=VOLUME/2;

    assign(spinor2,mnl->pf,VOLUME/2);
    /* apply B */
    for(j = 0; j < mnl->MDPolyDegree/2; j++){
      assign(spinor1, spinor2, VOLUME/2);
      Poly_Qtm_pm_min_cconst_nrm(spinor2,
			    spinor1,
			    mnl->MDPolyRoots[j]);
    }
    mnl->energy1 =  square_norm(spinor2, VOLUME/2,1);

    /* needed for return check */
/*     assign(mnl->pf,spinor2,VOLUME/2); */

    /* calculate evs */

    if (compute_evs != 0) {
      no_eigenvalues=10;
      eigenvalues(&no_eigenvalues, max_solver_iterations, eigenvalue_precision,
                  0/* compute minimal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);

      no_eigenvalues=1;
      eigenvalues(&no_eigenvalues, max_solver_iterations, eigenvalue_precision,
                  1/* compute maximal evs*/, 0/*dont write evecs*/, nstore, mnl->even_odd_flag);
    }


    
    g_mu = g_mu1;
    boundary(g_kappa);
    /* return the energy differnce */
    
    if(g_proc_id==0) {
      fprintf(stderr," Poly energy1     = %e \n" , mnl->energy1);
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
      diff = mnl->energy1 - mnl->energy0;
      fprintf(stderr," Poly energy diff = %e \n" , diff);
    }
    
    return (mnl->energy1 - mnl->energy0);
  } else {
    if(g_proc_id==0){
      fprintf(stderr,"Error: PHMC for light quarks not implementeted for non even/odd preconditioning\n");
    }

    g_mu = g_mu1;
    boundary(g_kappa);

    return NAN;
  }

}

void poly_heatbath(const int id){
  monomial * mnl = &monomial_list[id];
  g_mu = mnl->mu;
  boundary(mnl->kappa);
  int j;
  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;
  double sqdiff=-1;

  spinor* spinor1=g_spinor_field[2];
  spinor* spinor2=g_spinor_field[3];
  spinor* spinor3=g_spinor_field[4];

  actual_mnl=mnl;

  if(mnl->even_odd_flag) {
    poly_monomial_V=VOLUME/2;

/*     printf(stderr," testing whether P_n is kind of inverse .... \n"); */

/*     random_spinor_field(spinor1, VOLUME/2, mnl->rngrepro); */
/*     Poly_Ptm_pm_psi(spinor2,spinor1); */

/*     Qtm_pm_psi(spinor3,spinor2); */
/*     mul_r(spinor3, 1.0/ mnl->MDPolyLmax, spinor3, VOLUME/2); */


/*     diff(spinor2,spinor3,spinor1,VOLUME/2); */
/*     sqdiff=square_norm(spinor2,VOLUME/2,1); */
/*     fprintf(stderr,"  | ( P(D) * D - 1 ) \\phi_rnd |^2  = %e \n .... [done]\n", sqdiff ); */

    random_spinor_field(spinor1, VOLUME/2, mnl->rngrepro);
    mnl->energy0 = square_norm(spinor1, VOLUME/2, 1);

    if(g_proc_id==0) {
      fprintf(stderr," Poly energy0     = %e \n" , mnl->energy0);
    }



      /* calculate the phmc hamiltonian */
      Qtm_pm_psi(spinor2, spinor1);


/*       mnl->iter0 = cg_her(g_spinor_field[DUM_DERI+5], mnl->pf,  */
/* 			  mnl->maxiter, mnl->accprec, g_relative_precision_flag,  */
/* 			  VOLUME, Q_pm_psi, 0, 0); */


      /* solve (Q+)*(Q-)*P((Q+)*(Q-)) *x=y */
      /* TODO: define Qtm_pm_Ptm_pm_psi */
      cg_her(spinor1, spinor2,
             1000,mnl->accprec,g_relative_precision_flag,VOLUME/2, Poly_Qtm_pm_Ptm_pm_psi  ,0,1);

      /*  phi= Bdagger phi  */
      for(j = 0; j < (mnl->MDPolyDegree/2); j++){
	assign(spinor2, spinor1, VOLUME/2);
	Poly_Qtm_pm_min_cconst_nrm(spinor1,
			      spinor2,
			      mnl->MDPolyRoots[mnl->MDPolyDegree/2+j]);
      }
      assign(mnl->pf, spinor1, VOLUME/2);


/*     Qtm_plus_psi(mnl->pf, g_spinor_field[2]); */
/*     chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array, */
/* 			mnl->csg_N, &mnl->csg_n, VOLUME/2); */
/*     if(mnl->solver != CG) { */
/*       chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2,  */
/* 			  mnl->csg_N2, &mnl->csg_n2, VOLUME/2); */


    
  }
  else {
    if(g_proc_id==0)
      poly_monomial_c_msg("Error: NO NON-EVEN-ODD PRECONDITIONING IMPLEMENTED FOR PHMC");
  }

  g_mu = g_mu1;
  boundary(g_kappa);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called poly_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }

  return;
}


