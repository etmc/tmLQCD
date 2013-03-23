/***********************************************************************
 *
 * Copyright (C) 2006 Thomas Chiarappa
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
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_nd.h"
#include "chebyshev_polynomial_nd.h"
#include "phmc.h"
#include "solver/matrix_mult_typedef_nd.h"
#include "Ptilde_nd.h"


#define PI 3.141592653589793


double func_tilde(double u, double exponent){

  double ff=0.0;
  double d=0,ddd=0, sv, z, z2;
  int j;
  double res=0.0;

  z = (2.0*u - phmc_cheb_evmin - phmc_cheb_evmax)/(double)(phmc_cheb_evmax - phmc_cheb_evmin);
  z2 = 2.0*z;

  for(j=phmc_dop_n_cheby-1; j>=1; j--){
    sv = d;
    d = z2*d - ddd + phmc_dop_cheby_coef[j];
    ddd = sv;
  }

  res = z*d - ddd + 0.5*phmc_dop_cheby_coef[0];

  ff = (double)(res * sqrt(u));

  return(pow(ff,exponent));
}

void Ptilde_cheb_coefs(double aa, double bb, double dd[], int n, double exponent){
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
    f[k]=func_tilde(y*bma+bpa,exponent);
  }
  fac=2.0*inv_n;
  for (j=0;j<n;j++) {
    double sum=0.0;
    for (k=0;k<n;k++)
      sum += f[k]*cos(PI*j*(k+0.5)*inv_n);      
    dd[j]=fac*sum;
  }
  free(f);
}
#undef PI


/****************************************************************************  
 *
 * computation of the second poplynomial, Ptilde, on a vector
 *   by using the chebyshev approximation for the function ()^1/4
 * subtraction of low-lying eigenvalues is not yet implemented for this
 *
 **************************************************************************/

void Ptilde_ndpsi(spinor *R_s, spinor *R_c, double *dd, int n, 
		  spinor *S_s, spinor *S_c, matrix_mult_nd Qsq) {
  
  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4;
  
  spinor *svs_=NULL, *svs=NULL, *ds_=NULL, *ds=NULL, *dds_=NULL, *dds=NULL, 
    *auxs_=NULL, *auxs=NULL, *aux2s_=NULL, *aux2s=NULL, *aux3s_=NULL, 
    *aux3s=NULL;
  spinor *svc_=NULL, *svc=NULL, *dc_=NULL, *dc=NULL, *ddc_=NULL, 
    *ddc=NULL, *auxc_=NULL, *auxc=NULL, *aux2c_=NULL, *aux2c=NULL, 
    *aux3c_=NULL, *aux3c=NULL;
  
  
  svs_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  svs   = (spinor *)(((unsigned long int)(svs_)+ALIGN_BASE)&~ALIGN_BASE);
  ds_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  ds    = (spinor *)(((unsigned long int)(ds_)+ALIGN_BASE)&~ALIGN_BASE);
  dds_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  dds   = (spinor *)(((unsigned long int)(dds_)+ALIGN_BASE)&~ALIGN_BASE);
  auxs_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  auxs  = (spinor *)(((unsigned long int)(auxs_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2s_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux2s = (spinor *)(((unsigned long int)(aux2s_)+ALIGN_BASE)&~ALIGN_BASE);
  aux3s_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux3s = (spinor *)(((unsigned long int)(aux3s_)+ALIGN_BASE)&~ALIGN_BASE);
  svc_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  svc   = (spinor *)(((unsigned long int)(svc_)+ALIGN_BASE)&~ALIGN_BASE);
  dc_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  dc    = (spinor *)(((unsigned long int)(dc_)+ALIGN_BASE)&~ALIGN_BASE);
  ddc_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  ddc   = (spinor *)(((unsigned long int)(ddc_)+ALIGN_BASE)&~ALIGN_BASE);
  auxc_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  auxc  = (spinor *)(((unsigned long int)(auxc_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2c_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux2c = (spinor *)(((unsigned long int)(aux2c_)+ALIGN_BASE)&~ALIGN_BASE);
  aux3c_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux3c = (spinor *)(((unsigned long int)(aux3c_)+ALIGN_BASE)&~ALIGN_BASE);
  
  fact1=4/(phmc_cheb_evmax-phmc_cheb_evmin);
  fact2=-2*(phmc_cheb_evmax+phmc_cheb_evmin)/(phmc_cheb_evmax-phmc_cheb_evmin);
  
  zero_spinor_field(&ds[0],VOLUME/2);
  zero_spinor_field(&dds[0],VOLUME/2); 
  zero_spinor_field(&dc[0],VOLUME/2);
  zero_spinor_field(&ddc[0],VOLUME/2); 
  
  /*   sub_low_ev(&aux3[0], &S[0]);  */
  assign(&aux3s[0], &S_s[0],VOLUME/2);  
  assign(&aux3c[0], &S_c[0],VOLUME/2);  
  
  /*  Use the Clenshaw's recursion for the Chebysheff polynomial */
  for (j=n-1; j>=1; j--) {
    assign(&svs[0],&ds[0],VOLUME/2);
    assign(&svc[0],&dc[0],VOLUME/2); 
    
    /*
     * if ( (j%10) == 0 ) {
     *   sub_low_ev(&aux[0], &d[0]);
     * } else { */
    assign(&auxs[0], &ds[0], VOLUME/2);
    assign(&auxc[0], &dc[0], VOLUME/2);
    /*   } */


    Qsq(&R_s[0], &R_c[0], &auxs[0], &auxc[0]);

    temp1=-1.0;
    temp2=dd[j];
    assign_mul_add_mul_add_mul_add_mul_r(&ds[0] , &R_s[0], &dds[0], &aux3s[0], fact2, fact1, temp1, temp2,VOLUME/2);
    assign_mul_add_mul_add_mul_add_mul_r(&dc[0] , &R_c[0], &ddc[0], &aux3c[0], fact2, fact1, temp1, temp2,VOLUME/2);
    assign(&dds[0], &svs[0],VOLUME/2);
    assign(&ddc[0], &svc[0],VOLUME/2);
  }

  assign(&R_s[0], &ds[0],VOLUME/2);
  assign(&R_c[0], &dc[0],VOLUME/2);

  Qsq(&auxs[0], &auxc[0], &R_s[0], &R_c[0]);

  temp1=-1.0;
  temp2=dd[0]/2;
  temp3=fact1/2;
  temp4=fact2/2;
  assign_mul_add_mul_add_mul_add_mul_r(&auxs[0], &ds[0], &dds[0], &aux3s[0], temp3, temp4, temp1, temp2,VOLUME/2);
  assign_mul_add_mul_add_mul_add_mul_r(&auxc[0], &dc[0], &ddc[0], &aux3c[0], temp3, temp4, temp1, temp2,VOLUME/2);
  assign(&R_s[0], &auxs[0],VOLUME/2);
  assign(&R_c[0], &auxc[0],VOLUME/2);

  free(svs_);
  free(ds_);
  free(dds_);
  free(auxs_);
  free(aux2s_);
  free(aux3s_);
  free(svc_);
  free(dc_);
  free(ddc_);
  free(auxc_);
  free(aux2c_);
  free(aux3c_);
}

double chebtilde_eval(int M, double *dd, double s){

  double d=0,ddd=0, sv, z, z2, res;
  int j;

  z = (2.0*s - phmc_cheb_evmin - phmc_cheb_evmax)/(double)(phmc_cheb_evmax - phmc_cheb_evmin);
  z2 = 2.0*z;

  for(j=M-1; j>=1; j--){
    sv = d;
    d = z2*d - ddd + dd[j];
    ddd = sv;
  }

  res = z*d - ddd + 0.5*dd[0];

  return(res);
}

/**************************************************************************
 *
 * The externally accessible function is
 *
 *   void degree_of_Ptilde
 *     Computation of (QdaggerQ)^1/4
 *     by using the chebyshev approximation for the function ()^1/4  
 *
 * Author: Thomas Chiarappa <Thomas.Chiarappa@mib.infn.it> May 2006 
 *
 *****************************************************************************/



void degree_of_Ptilde(int * _degree, double ** coefs,
		      const double EVMin, const double EVMax,
		      const int sloppy_degree, const double acc, 
		      matrix_mult_nd Qsq, const int repro) {
  int i, j;
  double temp, temp2;
  int degree;
  double sum=0.0;

  spinor *ss=NULL, *ss_=NULL, *sc=NULL, *sc_=NULL;
  spinor *auxs=NULL, *auxs_=NULL, *auxc=NULL, *auxc_=NULL;
  spinor *aux2s=NULL, *aux2s_=NULL, *aux2c=NULL, *aux2c_=NULL;

  *coefs = calloc(phmc_max_ptilde_degree, sizeof(double)); 

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

  Ptilde_cheb_coefs(EVMin, EVMax, *coefs, phmc_max_ptilde_degree, -1.0); 

  if(g_proc_id == g_stdio_proc && g_debug_level > 0){
    printf("# NDPOLY Acceptance Polynomial: EVmin = %f  EVmax = %f\n", EVMin, EVMax);
    printf("# NDPOLY ACceptance Polynomial: desired accuracy is %e \n", acc);
    fflush(stdout);
  }

  degree = 2*sloppy_degree;

  for(i = 0; i < 100 ; i++) {
    if (degree > phmc_max_ptilde_degree) {
      fprintf(stderr, "Error: n_cheby=%d > phmc_max_ptilde_degree=%d in ptilde\n",
              degree, phmc_max_ptilde_degree);
      fprintf(stderr, "Increase n_chebymax\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(-5);
    }

    sum=0;
    for(j=degree; j<phmc_max_ptilde_degree; j++){ 
      sum += fabs(coefs[0][j]);
    }

    if((g_proc_id == g_stdio_proc) && (g_debug_level > 0)) {
      printf("# NDPOLY Acceptance Polynomial: Sum remaining | d_n | = %e for degree=%d\n", sum, degree);
      printf("# NDPOLY Acceptance Polynomial: coef[degree] = %e\n", (*coefs)[degree]);
    }
    if(sum < acc) { 
      break;
    }
    degree= (int)(degree*1.2);
  }

  if(g_debug_level > 2) {
    /* Ptilde P S P  Ptilde X - X */
    /* for random spinor X        */
    random_spinor_field_eo(ss, repro, RN_GAUSS);
    random_spinor_field_eo(sc, repro, RN_GAUSS);

    Ptilde_ndpsi(&auxs[0], &auxc[0], *coefs, degree, &ss[0], &sc[0], Qsq);
    Ptilde_ndpsi(&aux2s[0], &aux2c[0], phmc_dop_cheby_coef, phmc_dop_n_cheby, &auxs[0], &auxc[0], Qsq);
    Qsq(&auxs[0], &auxc[0], &aux2s[0], &aux2c[0]);
    Ptilde_ndpsi(&aux2s[0], &aux2c[0], phmc_dop_cheby_coef, phmc_dop_n_cheby, &auxs[0], &auxc[0], Qsq);
    Ptilde_ndpsi(&auxs[0], &auxc[0], *coefs, degree, &aux2s[0], &aux2c[0], Qsq);

    diff(&aux2s[0],&auxs[0], &ss[0], VOLUME/2);
    temp = square_norm(&aux2s[0], VOLUME/2, 1) / square_norm(&ss[0], VOLUME/2, 1) / 4.0;

    diff(&aux2c[0],&auxc[0], &sc[0], VOLUME/2);
    temp2 = square_norm(&aux2c[0], VOLUME/2, 1)/square_norm(&sc[0], VOLUME/2, 1) / 4.0;

    if(g_epsbar == 0){
      temp2 = 0.0;
    }
    /* || (Ptilde P S P Ptilde - 1)X ||^2 / || 2X ||^2 */
    if(g_proc_id == g_stdio_proc) {
      printf("# NDPOLY Acceptance Polynomial: relative squared accuracy in components:\n# UP=%e  DN=%e \n", temp, temp2);
    }

    temp = chebtilde_eval(degree, *coefs, EVMin);
    temp *= cheb_eval(phmc_dop_n_cheby, phmc_dop_cheby_coef, EVMin);
    temp *= EVMin;
    temp *= cheb_eval(phmc_dop_n_cheby, phmc_dop_cheby_coef, EVMin);
    temp *= chebtilde_eval(degree, *coefs, EVMin);
    temp = 0.5*fabs(temp - 1);
    if(g_proc_id == g_stdio_proc) {
      printf("# NDPOLY Acceptance Polynomial: Delta_IR at s=%f: | Ptilde P s_low P Ptilde - 1 |/2 = %e \n", EVMin, temp);
    }
  }
  if(g_proc_id == g_stdio_proc) {
    printf("# NDPOLY Acceptance Polynomial degree set to %d\n\n", degree);
  }

  *_degree = degree;
  free(ss_);
  free(auxs_);
  free(aux2s_);
  free(sc_);
  free(auxc_);
  free(aux2c_);
  return;
}
