/***********************************************************************
 * 
 * Copyright (C) 2003,2005,2006,2007,2008 Mauro Papinutto, Ines Wetzorke,
 *                                        Karl Jansen, Carsten Urbach
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
#include "chebyshev_polynomial.h"

#define PI 3.141592653589793

double cheb_evmin, cheb_evmax;

double func(double u, double exponent){
  return pow(u,exponent);
}

void chebyshev_polynomial(double aa, double bb, double c[], int n, double exponent){
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


/****************************************************************************  
 *
 * computation of (Q^dagger Q)^(1/4) on a vector
 *   by using the chebyshev approximation for the function ()^1/4
 * subtraction of low-lying eigenvalues is not yet implemented for this
 *
 **************************************************************************/

void QdaggerQ_power(spinor *R_s, spinor *R_c, double *c, int n, spinor *S_s, spinor *S_c)
{
  int j;
  double fact1, fact2, temp, temp1, temp2, temp3, temp4;
  spinor *svs_=NULL, *svs=NULL, *ds_=NULL, *ds=NULL, *dds_=NULL, *dds=NULL, 
 *auxs_=NULL, *auxs=NULL, *aux2s_=NULL, *aux2s=NULL, *aux3s_=NULL, *aux3s=NULL;
  spinor *svc_=NULL, *svc=NULL, *dc_=NULL, *dc=NULL, *ddc_=NULL, 
  *ddc=NULL, *auxc_=NULL, 
  *auxc=NULL, *aux2c_=NULL, *aux2c=NULL, *aux3c_=NULL, *aux3c=NULL;

  


#if ( defined SSE || defined SSE2 )
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
#else
   svs_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   svs = svs_;
   ds_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   ds = ds_;
   dds_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   dds = dds_;
   auxs_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   auxs = auxs_;
   aux2s_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   aux2s = aux2s_;
   aux3s_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   aux3s = aux3s_;
   svc_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   svc = svc_;
   dc_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   dc = dc_;
   ddc_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   ddc = ddc_;
   auxc_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   auxc = auxc_;
   aux2c_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   aux2c = aux2c_;
   aux3c_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
   aux3c = aux3c_;
#endif
 
   cheb_evmax=1.;


   fact1=4/(cheb_evmax-cheb_evmin);
   fact2=-2*(cheb_evmax+cheb_evmin)/(cheb_evmax-cheb_evmin);
   
   zero_spinor_field(&ds[0],VOLUME/2);
   zero_spinor_field(&dds[0],VOLUME/2); 
   zero_spinor_field(&dc[0],VOLUME/2);
   zero_spinor_field(&ddc[0],VOLUME/2); 
   
/*   sub_low_ev(&aux3[0], &S[0]);  */ 
    assign(&aux3s[0], &S_s[0],VOLUME/2);  
    assign(&aux3c[0], &S_c[0],VOLUME/2);  
   
     /* Use the Clenshaw's recursion for the 
	Chebysheff polynomial 
     */


     for (j=n-1; j>=1; j--) {
       assign(&svs[0],&ds[0],VOLUME/2); 
       assign(&svc[0],&dc[0],VOLUME/2); 
       
/*       if ( (j%10) == 0 ) {
  	 sub_low_ev(&aux[0], &d[0]);
         }
         else {    */
	 assign(&auxs[0], &ds[0], VOLUME/2);
	 assign(&auxc[0], &dc[0], VOLUME/2);
/*       }  */
       
       Qtm_dagger_ndpsi(&aux2s[0], &aux2c[0], &auxs[0], &auxc[0]);
       Qtm_ndpsi(&R_s[0], &R_c[0], &aux2s[0], &aux2c[0]);
       temp1=-1.0;
       temp2=c[j];
       assign_mul_add_mul_add_mul_add_mul_r(&ds[0] , &R_s[0], &dds[0], &aux3s[0], fact2, fact1, temp1, temp2,VOLUME/2);
       assign_mul_add_mul_add_mul_add_mul_r(&dc[0] , &R_c[0], &ddc[0], &aux3c[0], fact2, fact1, temp1, temp2,VOLUME/2);
       assign(&dds[0], &svs[0],VOLUME/2);
       assign(&ddc[0], &svc[0],VOLUME/2);
     } 
     
/*     sub_low_ev(&R[0],&d[0]);  */ 
     assign(&R_s[0], &ds[0],VOLUME/2);   
     assign(&R_c[0], &dc[0],VOLUME/2);  
     
     Qtm_dagger_ndpsi(&aux2s[0], &aux2c[0], &R_s[0], &R_c[0]);
     Qtm_ndpsi(&auxs[0], &auxc[0], &aux2s[0], &aux2c[0]);

     temp1=-1.0;
     temp2=c[0]/2;
     temp3=fact1/2;
     temp4=fact2/2;
     assign_mul_add_mul_add_mul_add_mul_r(&auxs[0], &ds[0], &dds[0], &aux3s[0], temp3, temp4, temp1, temp2,VOLUME/2);
     assign_mul_add_mul_add_mul_add_mul_r(&auxc[0], &dc[0], &ddc[0], &aux3c[0], temp3, temp4, temp1, temp2,VOLUME/2);
     assign(&R_s[0], &auxs[0],VOLUME/2);   
     assign(&R_c[0], &auxc[0],VOLUME/2);  
     
/*     addproj_q_invsqrt(&R[0], &S[0]); */
    
/*
#ifndef _SOLVER_OUTPUT
     if(g_proc_id == g_stdio_proc){
       printf("Order of Chebysheff approximation = %d\n",j); 
       fflush( stdout);};
#endif
*/
     
    
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
  


/**************************************************************************
 *
 * The externally accessible function is
 *
 *   void degree_of_polynomial(void)
 *     Computation of (QdaggerQ)^1/4
 *     by using the chebyshev approximation for the function ()^1/4  
 *
 *
 *******************************************************************************/

double stopeps=5.0e-16;

int dop_n_cheby=0;
double * dop_cheby_coef;

void degree_of_polynomial(const int repro){
  int i;
  double temp;
  static int ini=0;

  spinor *ss=NULL, *ss_=NULL, *auxs=NULL, *auxs_=NULL, 
         *aux2s=NULL, *aux2s_=NULL, *aux3s=NULL, *aux3s_=NULL;
  spinor *sc=NULL, *sc_=NULL, *auxc=NULL, *auxc_=NULL, *aux2c=NULL, 
         *aux2c_=NULL, *aux3c=NULL, *aux3c_=NULL;



  if(ini==0){
    dop_cheby_coef = calloc(N_CHEBYMAX,sizeof(double));
    ini=1;
  }




#if ( defined SSE || defined SSE2 || defined SSE3)
   ss_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   auxs_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   aux2s_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   aux3s_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   ss    = (spinor *)(((unsigned long int)(ss_)+ALIGN_BASE)&~ALIGN_BASE);
   auxs  = (spinor *)(((unsigned long int)(auxs_)+ALIGN_BASE)&~ALIGN_BASE);
   aux2s = (spinor *)(((unsigned long int)(aux2s_)+ALIGN_BASE)&~ALIGN_BASE);
   aux3s = (spinor *)(((unsigned long int)(aux3s_)+ALIGN_BASE)&~ALIGN_BASE);
   sc_   = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   auxc_ = calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   aux2c_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   aux3c_= calloc(VOLUMEPLUSRAND/2+1, sizeof(spinor));
   sc    = (spinor *)(((unsigned long int)(sc_)+ALIGN_BASE)&~ALIGN_BASE);
   auxc  = (spinor *)(((unsigned long int)(auxc_)+ALIGN_BASE)&~ALIGN_BASE);
   aux2c = (spinor *)(((unsigned long int)(aux2c_)+ALIGN_BASE)&~ALIGN_BASE);
   aux3c = (spinor *)(((unsigned long int)(aux3c_)+ALIGN_BASE)&~ALIGN_BASE);
#else
   ss   =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   auxs =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   aux2s=calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   aux3s=calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   sc   =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   auxc =calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   aux2c=calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
   aux3c=calloc(VOLUMEPLUSRAND/2, sizeof(spinor));
#endif

   chebyshev_polynomial(cheb_evmin, cheb_evmax, dop_cheby_coef, N_CHEBYMAX, 0.25);

   temp=1.0;
   random_spinor_field_eo(ss, repro, RN_GAUSS);
   random_spinor_field_eo(sc, repro, RN_GAUSS);
/*   assign(&sc[0], &ss[0],VOLUME/2);

  Qtm_pm_psi(&auxs[0], &ss[0]);
    temp=square_norm(&auxs[0],VOLUME/2, 1);
      printf("||auxs Carsten||=%e\n",temp);

  Qtm_dagger_ndpsi(&aux3s[0], &aux3c[0], &ss[0], &sc[0]);
  Qtm_ndpsi(&auxs[0], &auxc[0], &aux3s[0], &aux3c[0]);
    temp=square_norm(&auxs[0],VOLUME/2, 1);
      printf("||auxs own||=%e\n",temp);
    temp=square_norm(&auxc[0],VOLUME/2, 1);
      printf("||auxc own||=%e\n",temp); */


/*  if(g_proc_id == g_stdio_proc) {
    printf("\ndetermine the degree of the polynomial:\n");
    fflush(stdout);
  }  */

  dop_n_cheby=(int)5./sqrt(cheb_evmin);
  for(i = 0;i < 1 ; i++){
/*      printf("n_cheby=%d i=%d\n", dop_n_cheby, i); */
    
    if (dop_n_cheby >= N_CHEBYMAX) {
      if(g_proc_id == g_stdio_proc){
	printf("Error: n_cheby=%d > N_CHEBYMAX=%d\n",dop_n_cheby,N_CHEBYMAX);
	printf("Increase n_chebymax\n");
      }
/*      errorhandler(35,"degree_of_polynomial"); */ 
    }

    QdaggerQ_power(&aux3s[0], &aux3c[0], dop_cheby_coef, dop_n_cheby, &ss[0], &sc[0]);
    QdaggerQ_power(&auxs[0], &auxc[0], dop_cheby_coef, dop_n_cheby, &aux3s[0], &aux3c[0]);
    QdaggerQ_power(&aux3s[0], &aux3c[0], dop_cheby_coef, dop_n_cheby, &auxs[0], &auxc[0]);
    QdaggerQ_power(&auxs[0], &auxc[0], dop_cheby_coef, dop_n_cheby, &aux3s[0], &aux3c[0]);
/*    temp=square_norm(&auxs[0],VOLUME/2, 1);
    printf("||auxs||=%e\n",temp);
    temp=square_norm(&auxc[0],VOLUME/2, 1);
    printf("||auxc||=%e\n",temp); */


  Qtm_dagger_ndpsi(&aux2s[0], &aux2c[0], &ss[0], &sc[0]);
  Qtm_ndpsi(&aux3s[0], &aux3c[0], &aux2s[0], &aux2c[0]);

/*    temp=square_norm(&aux3s[0],VOLUME/2, 1);
      printf("||auxs_3||=%e\n",temp);
    temp=square_norm(&aux3c[0],VOLUME/2, 1);
      printf("||auxc_3||=%e\n",temp); */

    diff(&auxs[0],&auxs[0],&aux3s[0],VOLUME/2);
    temp=square_norm(&auxs[0],VOLUME/2)/square_norm(&aux3s[0],VOLUME/2, 1)/4.0;
    if(g_proc_id == g_stdio_proc) { 
      printf("difference=%e\n",temp);
    diff(&auxc[0],&auxc[0],&aux3c[0],VOLUME/2);
    temp=square_norm(&auxc[0],VOLUME/2)/square_norm(&aux3c[0],VOLUME/2, 1)/4.0;
      printf("difference=%e\n",temp);
   }  
    if(temp < stopeps ) break;
    dop_n_cheby*=1.05;
  }

   free(ss_);   
   free(auxs_); 
   free(aux2s_);
   free(aux3s_);
   free(sc_);   
   free(auxc_); 
   free(aux2c_);
   free(aux3c_);
}
