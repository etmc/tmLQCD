/***********************************************************************
 *
 * Copyright (C) 2011 Elena Garcia-Ramos
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
#include "su3.h"
#include "linalg_eo.h"
#include "chebyshev_polynomial_nd.h"
#include <io/eospinor.h>
#include "solver/solver.h"
#include "solver/jdher.h"
#include "solver/eigenvalues.h"
#include "X_psi.h"
#include "gamma.h"

double rnorm=-1;

/* |R>=rnorm^2 Q^2 |S> */
void norm_X_sqr_psi(spinor * const R, spinor * const S,
                    double const mstar);

/* |R>=rnorm Q|S> */
void norm_X_n_psi(spinor * const R, spinor * const S,
                  const int n, double const mstar);

/* Construct the sign function of the operator X  */
/* X/sqrt(X^2)   ,,   X = 1-(2M^2/(DdaggeraD+M^2))*/
void X_over_sqrt_X_sqr(spinor * const R, double * const c,
                       const int n, spinor * const S,
                       const double minev, double const mstar);


double * x_cheby_coef = NULL;
double epsilon=0.01; 
int x_n_cheby = 32;


void h_X_sqr_eta(spinor * const R1,spinor * const R2,spinor * const S, double const mstar){
  int i;
  double mode_n;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;
  
  /* Compute Chebyshev coefficients             */
  /* c[j] ,,  j=0..n, n=degree of the polynomial*/

  if(g_proc_id == 0) {
  printf("Degree of Polynomial set to %d\n", x_n_cheby);
  }

  if(n_cheby != x_n_cheby || rec_coefs) {
    if(x_cheby_coef != NULL) free(x_cheby_coef);
    x_cheby_coef = (double*)malloc(x_n_cheby*sizeof(double));
    chebyshev_coefs(epsilon, 1., x_cheby_coef, x_n_cheby, -0.5);//coefs for f(x)=x^(-0.5) // represents P(y)=1/sqrt(y) in paper "Chiral symmetry breaking an the Banks-Casher relation in lattice QCD with Wilson quarks" page 12.
    rec_coefs = 0;
    n_cheby = x_n_cheby;
  }

  if(g_proc_id == 0) {
  printf("mstar= %f \n",mstar);
  }

  /*Evaluate X_over_sqrt_X_sqr*/
  X_over_sqrt_X_sqr(R1, x_cheby_coef, x_n_cheby, S, epsilon, mstar);

  /* Construct h(x)=1/2-1/2 X/sqrt(X^2)                        */
  /* this routine makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
  assign_mul_add_mul_r(R1,S, 0.5, 0.5, VOLUME);

  
  /*we need h(X)^2|nu>*/
  X_over_sqrt_X_sqr(R2, x_cheby_coef, x_n_cheby, R1, epsilon, mstar);
  assign_mul_add_mul_r(R2,R1,0.5, 0.5, VOLUME);
  
  return;
}

void h_X_eta(spinor * const R,spinor * const S, double const mstar){
  int i;
  double mode_n;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  /* Compute Chebyshev coefficients             */
  /* c[j] ,,  j=0..n, n=degree of the polynomial*/

  if(g_proc_id == 0) {
  printf("Degree of Polynomial set to %d\n", x_n_cheby);
  }

  if(n_cheby != x_n_cheby || rec_coefs) {
    if(x_cheby_coef != NULL) free(x_cheby_coef);
    x_cheby_coef = (double*)malloc(x_n_cheby*sizeof(double));
    chebyshev_coefs(epsilon, 1., x_cheby_coef, x_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = x_n_cheby;
  }

  /*Evaluate X_over_sqrt_X_sqr*/
  X_over_sqrt_X_sqr(R, x_cheby_coef, x_n_cheby, S, epsilon, mstar);

  /* Construct h(x)=1/2-1/2 X/sqrt(X^2)                        */
  /* this routine makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
  assign_mul_add_mul_r(R,S, 0.5, 0.5, VOLUME);

  return;
}


void h_X_4_eta(spinor * const R1,spinor * const R2,spinor * const S, double const mstar){
  int i;
  double mode_n;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  /* Compute Chebyshev coefficients             */
  /* c[j] ,,  j=0..n, n=degree of the polynomial*/
  if(g_proc_id == 0) {
  printf("Degree of Polynomial set to %d\n", x_n_cheby);
  }

  if(n_cheby != x_n_cheby || rec_coefs) {
    if(x_cheby_coef != NULL) free(x_cheby_coef);
    x_cheby_coef = (double*)malloc(x_n_cheby*sizeof(double));
    chebyshev_coefs(epsilon, 1., x_cheby_coef, x_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = x_n_cheby;
  }
  s_ = calloc(3*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(3, sizeof(spinor*));

  for(i = 0; i < 3; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
  }

  printf("mstar= %f \n",mstar);
  /* Evaluate X_over_sqrt_X_sqr */
  X_over_sqrt_X_sqr(s[0], x_cheby_coef, x_n_cheby, S, epsilon, mstar);

  /* Construct h(x)=1/2-1/2 X/sqrt(X^2)                        */
  /* this routine makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
  assign_mul_add_mul_r(s[0],S, 0.5, 0.5, VOLUME);

  X_over_sqrt_X_sqr(R1, x_cheby_coef, x_n_cheby, s[0], epsilon, mstar);
  assign_mul_add_mul_r(R1,s[0],0.5, 0.5, VOLUME);

  X_over_sqrt_X_sqr(s[2], x_cheby_coef, x_n_cheby, R1, epsilon, mstar);
  assign_mul_add_mul_r(s[2],R1,0.5, 0.5, VOLUME);

  /*we need h(X)^2|nu>*/
  X_over_sqrt_X_sqr(R2, x_cheby_coef, x_n_cheby, s[2], epsilon, mstar);
  assign_mul_add_mul_r(R2,s[2],0.5, 0.5, VOLUME);

  free(s);
  free(s_);

  return;
}



void norm_X_sqr_psi(spinor * const R, spinor * const S, double const mstar) {

  spinor *aux_,*aux;
#if ( defined SSE || defined SSE2 || defined SSE3 )
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif

  /* Here is where we have to include our operator which in this case is
     X = 1 - (2M^2)/(D_m^dagger*D_m + mu^2 + M^2)  */
  
  if(1)
  {
    X_psi(aux, S, mstar);
    X_psi(R, aux, mstar);
  }
  else
  {
    printf("using X_psiSquare.\n");
    X_psiSquare(R, S, mstar);
  }
  mul_r(R, rnorm*rnorm, R, VOLUME);


  free(aux_);
  return;
}


void norm_X_n_psi(spinor * const R, spinor * const S, 
		  const int n, double const mstar) { 

  int i;
  double npar = 1.;
  spinor *aux_,*aux;
#if (defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif
  assign(aux, S, VOLUME);
  
  for(i=0; i < n; i++){
  /* Here is where we have to include our operator which in this case is
     X = 1 - (2M^2)/(D_m^dagger*D_m + M^2)  */
    X_psi(R, aux, mstar);
    npar *= rnorm;
  }
  mul_r(R, npar, R, VOLUME);

  free(aux_);
  return;
}

void X_over_sqrt_X_sqr(spinor * const R, double * const c, 
		       const int n, spinor * const S, const double minev, double const mstar) {
//x/sqrt(x*x) <=> normalisation <= reasoned by Clenshaw recurrence: maps X to [-1,1] 
 
  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4, maxev;
  spinor *sv_, *sv, *d_, *d, *dd_, *dd, *aux_, *aux, *aux3_, *aux3;// *_ holds the adress of the sse-unaligned memory block
  
#if ( defined SSE || defined SSE2 || defined SSE3)
  sv_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  sv   = (spinor *)(((unsigned long int)(sv_)+ALIGN_BASE)&~ALIGN_BASE);
  d_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  d    = (spinor *)(((unsigned long int)(d_)+ALIGN_BASE)&~ALIGN_BASE);
  dd_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  dd   = (spinor *)(((unsigned long int)(dd_)+ALIGN_BASE)&~ALIGN_BASE);
  aux_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux  = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
  aux3_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux3 = (spinor *)(((unsigned long int)(aux3_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  sv_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  sv = sv_;
  d_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  d = d_;
  dd_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  dd = dd_;
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
  aux3_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux3 = aux3_;
#endif

  /*EVALUATE THE APPROXIMATION USING THE CLENSHAW'S RECURRENCE FORMULA*/

  maxev=1.0;

  /*interval = [minev,maxev] = [epsilon,1]*/  
  fact1=4/(maxev-minev);
  fact2=-2*(maxev+minev)/(maxev-minev);
  /* d=0 , dd=0 */
  zero_spinor_field(d, VOLUME);
  zero_spinor_field(dd, VOLUME); 

  
  /*input S = aux3*/
  if(0) assign_sub_lowest_eigenvalues(aux3, S, no_eigenvalues-1, VOLUME);
  else assign(aux3, S, VOLUME);


  /*starting the loop*/
  if(1) {
    for (j = n-1; j >= 1; j--) {

      /*sv=d  = d_j+1*/
      assign(sv, d, VOLUME); 

      /*aux= our random field S =0(j=n-1)*/
      assign(aux, d, VOLUME);

      if(j == n-1){
	assign(R, aux, VOLUME);//=0
      }
      else{
      /*|R>=rnorm^2 X^2|aux> -> since aux=d -> |R>=rnorm^2 Q^2|d>*/
        norm_X_sqr_psi(R, aux, mstar);//WARNING: - maybe we have to pass this point only when j=n-2, because R is not manipulated  in the loop body.
                                      //         - seems to setup d_n-1=0
      }
      temp1=-1.0;
      temp2=c[j]; /*Chebyshev coefficients*/

      /* d = d*fact2 + R*fact1 + dd*temp1 + aux3*temp2 
         d = -2*(maxev+minev)/(maxev-minev)*d + 4/(maxev-minev)*R 
	      -1*dd + c[j]*aux3                                           */
      /* y = (2*x-a-b)/(b-a)   ,   y2=2*y 
         d = y2*d - dd + c[j] = -2*(a+b)*d/(b-a) + 4*x*d/(b-a) -dd + c[j] */
      assign_mul_add_mul_add_mul_add_mul_r(d, R, dd, aux3, fact2, fact1, temp1, temp2, VOLUME);// =d_j+1
      /* dd = sv */
      assign(dd, sv, VOLUME);// = d_j+2
    } 
    
    /* R = d */
    if(0) assign_sub_lowest_eigenvalues(R, d, no_eigenvalues-1, VOLUME);
    else assign(R, d, VOLUME);
    
    /*|aux>=rnorm^2 Q^2|R> */
    norm_X_sqr_psi(aux, R, mstar);
    temp1=-1.0;
    temp2=c[0]/2.;
    temp3=fact1/2.;
    temp4=fact2/2.;
    
    /* aux = aux*temp3 + d*temp4 + dd*temp1 + aux3*temp2
       aux = 2/(maxev-minev)*aux + -(maxev+minev)/(maxev-minev)d 
             -1*dd + 0.5*c[j]*aux3                                */
    /* P(X^2)|_x = y*d -dd + 0.5*c[0]                             */
    assign_mul_add_mul_add_mul_add_mul_r(aux, d, dd, aux3, temp3, temp4, temp1, temp2, VOLUME);
    /* ONCE WE HAVE THE EVALUATION OF P(X^2) = 1/SQRT(X^2)
        WE CONSTRUCT -X/SQRT(X^2) -->  -X*P(X^2) */
    norm_X_n_psi(R, aux, 1, mstar);
  }
  
  free(sv_);
  free(d_);
  free(dd_);
  free(aux_);
  free(aux3_);
  return;
}


void Check_Approximation(double const mstar, const int repro) {

  if(g_proc_id == 0) {
  printf("Checking the approximation of X/sqrt(X^2) in the mode number: \n");
  }

  int i;
  double res = 0;
  spinor **s, *s_;
  spinor *Sin = NULL;
  spinor *Sin_ = NULL;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  printf("epsilon= %f \n", epsilon);
  printf("M*^2= %f \n", mstar);
  printf("x_n_cheby= %d \n", x_n_cheby);
  if(n_cheby != x_n_cheby || rec_coefs) {
    if(x_cheby_coef != NULL) free(x_cheby_coef);
    x_cheby_coef = (double*)malloc(x_n_cheby*sizeof(double));
    chebyshev_coefs(epsilon, 1., x_cheby_coef, x_n_cheby, -0.5);
    rec_coefs = 0;
    n_cheby = x_n_cheby;
  }
  
#if (defined SSE3 || defined SSE2 || defined SSE)
  Sin_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  Sin    = (spinor *)(((unsigned long int)(Sin_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  Sin   =calloc(VOLUMEPLUSRAND, sizeof(spinor));
#endif

  random_spinor_field_lexic(Sin, repro, RN_GAUSS);

  s_ = calloc(4*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(4, sizeof(spinor*));

  for(i = 0; i < 4; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
  }

  X_over_sqrt_X_sqr(s[0], x_cheby_coef, x_n_cheby, Sin, epsilon, mstar);
  
  diff(s[2], Sin, s[0], VOLUME);
  diff(s[2], Sin, s[0], VOLUME);
  
  X_over_sqrt_X_sqr(s[1], x_cheby_coef, x_n_cheby, s[0], epsilon, mstar);
  
  diff(s[3], s[1], Sin, VOLUME);
  res = square_norm(s[3],VOLUME,0);

  if(g_proc_id == 0) {
  printf("\n");
  printf("Deviation from the real value : \n");
  printf("||X^2/sqrt(X^2)|psi> - |nu>||^2 = %1.4e \n",res);
  printf("\n");
  }
  
#if (defined SSE3 || defined SSE2 || defined SSE)
  free(Sin_);
#else
  free(Sin);
#endif
  free(s);
  free(s_);
  return;
}

