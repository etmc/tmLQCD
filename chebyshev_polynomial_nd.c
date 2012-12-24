/***********************************************************************
 *
 * Copyright (C) 2006,2007,2008 Thomas Chiarappa, Carsten Urbach
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
#include "linalg_eo.h"
#include "start.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_nd.h"
#include "phmc.h"
#include "Ptilde_nd.h"
#include "chebyshev_polynomial_nd.h"



#define PI 3.141592653589793

double func(double u, double exponent){
  return pow(u,exponent);
}


void chebyshev_coefs(double aa, double bb, double c[], int n, double exponent){
  int k,j;
  double fac,bpa,bma,*f;
  double inv_n;


  inv_n=1./(double)n;
  f=calloc(n,sizeof(double));/*vector(0,n-1);*/
  fflush(stdout);
  bma=0.5*(bb-aa);
  bpa=0.5*(bb+aa);
  for (k=0;k<n;k++) {
    double y=cos(PI*(k+0.5)*inv_n);
    f[k]=func(y*bma+bpa,exponent);
  }
  fac=2.0*inv_n;
  for (j=0;j<n;j++) {
    double sum=0.0;
    for (k=0;k<n;k++)
      sum += f[k]*cos(PI*j*(k+0.5)*inv_n);
    c[j]=fac*sum;
  }
  free(f);


}
#undef PI


double cheb_eval(int M, double *c, double s){

  double d=0,dd=0, sv, z, z2, res;
  int j;

  z = (2.0*s - phmc_cheb_evmin - phmc_cheb_evmax)/(double)(phmc_cheb_evmax - phmc_cheb_evmin);
  z2 = 2.0*z;

  for(j=M-1; j>=1; j--){
    sv = d;
    d = z2*d - dd + c[j];
    dd = sv;
    }

  res = z*d - dd + 0.5*c[0];

  return(res);  
}

/**************************************************************************
 *
 * The externally accessible function is
 *
 *   void degree_of_polynomial_nd(void)
 *     Computation of (QdaggerQ)^1/4
 *     by using the chebyshev approximation for the function ()^1/4  
 *
 *
 *****************************************************************************/


void degree_of_polynomial_nd(int * _degree_of_p, double ** coefs,
			     const double EVMin, const double EVMax,
			     matrix_mult_nd Qsq, const int repro) { 
  double temp, temp2;
  int degree_of_p = *_degree_of_p + 1;

  spinor *ss=NULL, *ss_=NULL, *sc=NULL, *sc_=NULL;
  spinor *auxs=NULL, *auxs_=NULL, *auxc=NULL, *auxc_=NULL;
  spinor *aux2s=NULL, *aux2s_=NULL, *aux2c=NULL, *aux2c_=NULL;

  *coefs = calloc(degree_of_p, sizeof(double));

  ss_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  auxs_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  aux2s_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  sc_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  auxc_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  aux2c_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
  
  ss    = (spinor *)(((unsigned long int)(ss_)+ALIGN_BASE)&~ALIGN_BASE);
  auxs  = (spinor *)(((unsigned long int)(auxs_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2s = (spinor *)(((unsigned long int)(aux2s_)+ALIGN_BASE)&~ALIGN_BASE);
  sc    = (spinor *)(((unsigned long int)(sc_)+ALIGN_BASE)&~ALIGN_BASE);
  auxc  = (spinor *)(((unsigned long int)(auxc_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2c = (spinor *)(((unsigned long int)(aux2c_)+ALIGN_BASE)&~ALIGN_BASE);
  
  chebyshev_coefs(EVMin, EVMax, *coefs, degree_of_p, -0.5);

  random_spinor_field_eo(ss, repro, RN_GAUSS);
  random_spinor_field_eo(sc, repro, RN_GAUSS);

  if((g_proc_id == g_stdio_proc) && (g_debug_level > 0)){
    printf("# NDPOLY MD Polynomial: EVmin = %e  EVmax = %e  \n", EVMin, EVMax);
    printf("# NDPOLY MD Polynomial: the degree was set to: %d\n", degree_of_p);
    fflush(stdout);
  }

  if(g_debug_level > 1) {
    /* Here we check the accuracy */
    Ptilde_ndpsi(&auxs[0], &auxc[0], *coefs, degree_of_p, &ss[0], &sc[0], Qsq);
    Qsq(&aux2s[0], &aux2c[0], &auxs[0], &auxc[0]);
    Ptilde_ndpsi(&auxs[0], &auxc[0], *coefs, degree_of_p, &aux2s[0], &aux2c[0], Qsq);
    
    diff(&aux2s[0],&auxs[0],&ss[0],VOLUME/2);
    temp=square_norm(&aux2s[0],VOLUME/2, 1)/square_norm(&ss[0],VOLUME/2, 1)/4.0;
    
    diff(&aux2c[0],&auxc[0],&sc[0],VOLUME/2);
    temp2 = square_norm(&aux2c[0],VOLUME/2, 1)/square_norm(&sc[0],VOLUME/2, 1)/4.0;
    
    if(g_epsbar == 0.){ 
      temp2 = 0.0;
    }
    
    if(g_proc_id == g_stdio_proc){
      /* this is || (P S P - 1)X ||^2 /|| 2X ||^2 */
      /* where X is a random spinor field         */
      printf("# NDPOLY MD Polynomial: relative squared accuracy in components:\n# UP=%e  DN=%e \n", temp, temp2);
      fflush(stdout);
    }

    temp = cheb_eval(degree_of_p, *coefs, EVMin);
    temp *= EVMin;
    temp *= cheb_eval(degree_of_p, *coefs, EVMin);
    temp = 0.5*fabs(temp - 1);
    if(g_proc_id == g_stdio_proc) {
      printf("# NDPOLY MD Polynomial: Delta_IR at s=%f:    | P s_low P - 1 |/2 = %e \n", EVMin, temp);
    }
  }
  /* RECALL THAT WE NEED AN EVEN DEGREE !!!! */
  *_degree_of_p = degree_of_p;

  free(ss_);   
  free(auxs_); 
  free(aux2s_);
  free(sc_);   
  free(auxc_); 
  free(aux2c_);
  return;
}
