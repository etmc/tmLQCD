/* $Id$ */

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
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "linalg_eo.h"
#include "start.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "eigenvalues.h"
#include "derivative_psf.h"
#include "get_rectangle_staples.h"
#include "derivative_psf.h"
#include "gamma.h"
#include "get_staples.h"
#include "chebyshev_polynomial.h"
#include "Nondegenerate_Matrix.h"
#include "invert_eo.h"

#include "observables.h"
#include "hybrid_update.h"

#include "hybrid_nondegenerate_update.h"



/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

/* input is the pseudo-fermion field */
void deri_nondegenerate() {

  int i, mu, j, k;


  /* Recall:  The GAMMA_5 left of  delta M_eo  is done in  deriv_Sb !!! */


  /* Re-initialize df0 */
  for(i=0;i<(VOLUME+RAND);i++){ 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }


  /* Here comes the definitions for the chi_j fields */
  /* from  j=0  (chi_0 = phi)  .....  to j = n-1 */
  
  for(k=1; k<(dop_n_cheby-1); k++){

    L_POLY_MIN_CCONST(g_chi_up_spinor_field[k], g_chi_dn_spinor_field[k], g_chi_up_spinor_field[k-1], g_chi_dn_spinor_field[k-1], roo[k-1]);
  }


  /* Here comes the remaining fields  chi_k ; k=n,...,2n-1  */
  /*They are evaluated step-by-step overwriting the same field (dop_n_cheby)*/

  assign(g_chi_up_spinor_field[dop_n_cheby], g_chi_up_spinor_field[dop_n_cheby-2], VOLUME/2);
  assign(g_chi_dn_spinor_field[dop_n_cheby], g_chi_dn_spinor_field[dop_n_cheby-2], VOLUME/2);

  for(j=(dop_n_cheby-1); j>=1; j--){

    assign(g_chi_up_spinor_field[dop_n_cheby-1], g_chi_up_spinor_field[dop_n_cheby], VOLUME/2);
    assign(g_chi_dn_spinor_field[dop_n_cheby-1], g_chi_dn_spinor_field[dop_n_cheby], VOLUME/2);

    L_POLY_MIN_CCONST(g_chi_up_spinor_field[dop_n_cheby], g_chi_dn_spinor_field[dop_n_cheby], g_chi_up_spinor_field[dop_n_cheby-1], g_chi_dn_spinor_field[dop_n_cheby-1], roo[2*dop_n_cheby-3-j]);


    assign(g_spinor_field[DUM_DERI+4], g_chi_up_spinor_field[dop_n_cheby], VOLUME/2);
    assign(g_spinor_field[DUM_DERI+5], g_chi_dn_spinor_field[dop_n_cheby], VOLUME/2);
  
    assign(g_spinor_field[DUM_DERI+2], g_chi_up_spinor_field[j-1], VOLUME/2);
    assign(g_spinor_field[DUM_DERI+3], g_chi_dn_spinor_field[j-1], VOLUME/2);

      
    /* Get the even parts of the  (j-1)th  chi_spinors */
    H_eo_ND(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3], EO);

    /* \delta M_eo sandwitched by  chi[j-1]_e^\dagger  and  chi[2N-j]_o */
    deriv_Sb(EO, DUM_DERI, DUM_DERI+4);      /* UP */
    deriv_Sb(EO, DUM_DERI+1, DUM_DERI+5);    /* DN */


    /* Get the even parts of the  (2N-j)-th  chi_spinors */
    H_eo_ND(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+4], g_spinor_field[DUM_DERI+5], EO);

    /* \delta M_oe sandwitched by  chi[j-1]_o^\dagger  and  chi[2N-j]_e */
    deriv_Sb(OE, DUM_DERI+2, DUM_DERI);
    deriv_Sb(OE, DUM_DERI+3, DUM_DERI+1);
  }

  /*
    Normalisation by the largest  EW  is done in fermion_momenta_ND 
  */ 
}



void fermion_momenta_ND(double step) {
  int i,mu;
  double tmp;
  su3adj *xm,*deriv;

  deri_nondegenerate(); 

#ifdef MPI
  xchange_deri();
#endif
  for(i = 0; i < VOLUME; i++){
    for(mu=0;mu<4;mu++){
      xm=&moment[i][mu];
      
      _assign_const_times_mom(df0[i][mu], Cpol*invmaxev,df0[i][mu]);
      
      deriv=&df0[i][mu];
      /* Factor 2 around? */
      tmp = -2.*step;
      _minus_const_times_mom(*xm,tmp,*deriv); 
    }
   }

}


void leap_frog_ND(double step, int m, int nsmall) {
  int i,j;
  double smallstep;

  /* initialize the counter for the inverter */
  count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
  /* adjust the step-size to standard convention */
  step*=0.7071067811865;
  smallstep=step/nsmall;


#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta_ND(0.5*step);
  gauge_momenta(0.5*smallstep);
  for(i=1;i<m;i++){
    for(j=0;j<nsmall;j++){
      update_gauge(smallstep); 
      gauge_momenta(smallstep);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    fermion_momenta_ND(step);
  }
  for(j=1;j<nsmall;j++){
    update_gauge(smallstep); 
    gauge_momenta(smallstep);
  }
  update_gauge(smallstep); 
  gauge_momenta(0.5*smallstep);
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta_ND(0.5*step);

}


