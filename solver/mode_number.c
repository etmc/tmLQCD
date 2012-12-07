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
int x_n_cheby ;
double prec = 1.e-3;

void mode_number(spinor * const S, double const mstar){
  int i;
  double mode_n;
  spinor **s, *s_;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  //  x_n_cheby = (int)(-log(1.e-12)/(2*sqrt(epsilon)));
  x_n_cheby = (int)(-log(prec)/(2*sqrt(epsilon)));

  /*  if(g_proc_id == 0){
    printf("using precision %1.1e we get a degree n = %d \n",prec,x_n_cheby);
    }*/


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
  s_ = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(2, sizeof(spinor*));
  
  for(i = 0; i < 2; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
  }

  if(g_proc_id == 0) {
    printf("mstar= %f \n",mstar);
  }

  /*Evaluate X_over_sqrt_X_sqr*/
  X_over_sqrt_X_sqr(s[0], x_cheby_coef, x_n_cheby, S, epsilon, mstar);

  /* Construct h(x)=1/2-1/2 X/sqrt(X^2)                        */
  /* this routine makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
  assign_mul_add_mul_r(s[0],S, 0.5, 0.5, VOLUME);
  
  /*we need h(X)^2|nu>*/
  X_over_sqrt_X_sqr(s[1], x_cheby_coef, x_n_cheby, s[0], epsilon, mstar);
  assign_mul_add_mul_r(s[1],s[0],0.5, 0.5, VOLUME);
  
  /* Calculate the square norm ||h(X)^2 nu||^2 */
  mode_n = square_norm(s[1], VOLUME, 1); 

  if(g_proc_id == 0) {
    printf("The Value of the Mode Number is %f \n", mode_n);
  }


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
  
  X_psi(aux, S, mstar);
  X_psi(R, aux, mstar);
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
  
  int j;
  double fact1, fact2, temp1, temp2, temp3, temp4, maxev;
  spinor *sv_, *sv, *d_, *d, *dd_, *dd, *aux_, *aux, *aux3_, *aux3;
  //  double ap_eps_sq = 0.;
  
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

      /*sv=d*/
      assign(sv, d, VOLUME); 

      /*aux= our random field S*/
      assign(aux, d, VOLUME);

      if(j == n-1){
	assign(R, aux, VOLUME);
      }
      else{
      /*|R>=rnorm^2 X^2|aux> -> since aux=d -> |R>=rnorm^2 Q^2|d>*/
      norm_X_sqr_psi(R, aux, mstar);
      }
      temp1=-1.0;
      temp2=c[j]; /*Chebyshev coefficients*/

      /* d = d*fact2 + R*fact1 + dd*temp1 + aux*temp2 
         d = -2*(maxev+minev)/(maxev-minev)*d + 4/(maxev-minev)*R 
	      -1*dd + c[j]*aux3                                           */
      /* y = (2*x-a-b)/(b-a)   ,   y2=2*y 
         d = y2*d - dd + c[j] = -2*(a+b)*d/(b-a) + 4*x*d/(b-a) -dd + c[j] */
      assign_mul_add_mul_add_mul_add_mul_r(d, R, dd, aux3, fact2, fact1, temp1, temp2, VOLUME);
      /* dd = sv */
      assign(dd, sv, VOLUME);
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
    //, *Sin_ = NULL;
  static int n_cheby = 0;
  static int rec_coefs = 1;

  //  x_n_cheby = (int)(-log(1.e-12)/(2*sqrt(epsilon)));
  x_n_cheby = (int)(-log(prec)/(2*sqrt(epsilon)));

  if(g_proc_id == 0) {
  printf("epsilon= %f \n", epsilon);
  printf("M*^2= %f \n", mstar);
  printf("x_n_cheby= %d \n", x_n_cheby);
  }

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

  free(s);
  free(s_);
  return;
}

