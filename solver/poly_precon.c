/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "start.h"
#include "linalg_eo.h"
#include "tm_operators.h"
#include "chebyshev_polynomial.h"


#define PI 3.141592653589793

double f_pre(double u){
/*   return(1./(1.+(g_mu*g_mu*(1-1.5*1.5))/(u)));  */
/*   return pow(u,exponent); */
   return(1./u);
}

void get_c(double aa, double bb, double c[], int n){
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
    f[k]=f_pre(y*bma+bpa);
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


void poly_precon(spinor * const R, spinor * const S, const double prec, const int n) {
  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4, invmaxev = 1./4., maxev=4., tnorm, minev=g_mu*g_mu, auxnorm;
  static spinor *sv_, *sv, *d_, *d, *dd_, *dd, *aux_, *aux, *aux3_, *aux3;
  static int initp = 0;
  static double * c;
  const int N = VOLUME;


  
  maxev = 4.0;
  invmaxev = 1./maxev;
  minev = 0.1;
/*   minev = 1.5*1.5*g_mu*g_mu; */

  if(initp == 0) {
    c = (double*)calloc(1000, sizeof(double));
#if (defined SSE || defined SSE2 || defined SSE3)
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
    sv_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    sv   = sv_;
    d_   = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    d    = d_;
    dd_  = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    dd   = dd_;
    aux_ = calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    aux  = aux_;
    aux3_= calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
    aux3 = aux3_;
#endif
    get_c(minev, maxev, c, 100);
    initp = 1;
  }


  fact1 = 4. / (maxev - minev);
  fact2 = -2 * (maxev + minev) / (maxev - minev);
   
  zero_spinor_field(&d[0], N);
  zero_spinor_field(&dd[0], N); 
  assign(&aux3[0], &S[0], N); 
/*   gamma5(&aux3[0], &S[0], N); */

  /* Use the adaptive precision version using the forward recursion 
     for the Chebysheff polynomial 
  */

  /* d = T_0(Q^2) */
  assign(&d[0], &aux3[0], N);
  /* dd = T_1(Q^2) */
  Q_pm_psi(&dd[0], &d[0]);
/*   mul_r(dd, invmaxev, dd, N); */
  /*    norm_Q_sqr_psi(&dd[0], &d[0], g_m_D_psi, rnorm); */
  temp3 = fact1/2;
  temp4 = fact2/2;  
  assign_mul_add_mul_r(&dd[0], &d[0], temp3, temp4, N);
  /* r = c_1 T_1(Q^2) + 1/2 c_0 */
  temp1 = c[1];
  temp2 = c[0]/2;
  mul_add_mul_r(&R[0], &dd[0], &d[0], temp1, temp2, N);
     
  temp1 = -1.0;
  for (j=2; j<=n-1; j++) {
    /* aux = T_j(Q^2) = 2 Q^2 T_{j-1}(Q^2) - T_{j-2}(Q^2) */
    Q_pm_psi(&aux[0], &dd[0]);
/*     mul_r(aux, invmaxev, aux, N); */
    /*        norm_Q_sqr_psi(&aux[0], &dd[0], g_m_D_psi, rnorm); */
    assign_mul_add_mul_add_mul_r(&aux[0],&dd[0],&d[0],fact1,fact2,temp1, N);
    /* r = r + c_j T_j(Q^2) */
    temp2=c[j];
    assign_add_mul_r(&R[0],&aux[0],temp2, N);
    /* The stoppping criterio tnorm = |T_j(Q^2)| */
    tnorm = square_norm(aux, N);
    tnorm *= (temp2*temp2);
     
    
    auxnorm = square_norm(R, N);
    if(g_proc_id == g_stdio_proc) {
      printf("j= %d\t|c T|^2= %g\t%g\t c_j= %g\t|r|^2= %g\n",j,tnorm,prec, temp2,auxnorm); fflush( stdout);
      fflush(stdout);
    }
         
    if(tnorm < prec) break;
    /* d = T_{j-1}(Q^2) */
    assign(&d[0], &dd[0], N);
    /* dd = T_{j}(Q^2) */
    assign(&dd[0], &aux[0], N);
  }
  if(g_proc_id == g_stdio_proc) {
    printf("Order of Chebysheff approximation = %d\n",j); 
    fflush( stdout);
  }
   

  /* r = Q r */

/*   assign(aux, R, N); */
/*   Q_minus_psi(R, aux); */

  return;
}


