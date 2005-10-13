/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "eigenvalues.h"
#include "derivative_psf.h"
#include "gamma.h"
#include "chebyshev_polynomial.h"
#include "Nondegenerate_Matrix.h"
#include "invert_eo.h"
#include "hybrid_nondegenerate_update.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;


/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

/* input is the pseudo-fermion field */
void deri_nondegenerate() {

  int i,mu;

  for(i=0;i<(VOLUME+RAND);i++){ 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }

  /* This replaces the solver     */
  /* Apply the invers square root */
  /* to phi_o -> X_o              */
  QdaggerQ_power(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1], 
		 dop_cheby_coef, dop_n_cheby, 
		 spinor_field[STRANGE], spinor_field[CHARM]);
  

  /* Construct X_e and X_e' */
  /* First X_e              */
  /* normalization corrects?*/
  QNon_degenerate_eo(spinor_field[DUM_DERI+2], spinor_field[DUM_DERI+3], 
		     spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]); 
  mul_r(spinor_field[DUM_DERI+2], 1/20., spinor_field[DUM_DERI+2], VOLUME/2); 
  mul_r(spinor_field[DUM_DERI+3], 1/20., spinor_field[DUM_DERI+3], VOLUME/2); 

  /* X_e'                   */
  /* normalization corrects?*/
  QNon_degenerate_eo_dagger(spinor_field[DUM_DERI+4], spinor_field[DUM_DERI+5], 
		     spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]); 
  mul_r(spinor_field[DUM_DERI+4], 1/20., spinor_field[DUM_DERI+4], VOLUME/2); 
  mul_r(spinor_field[DUM_DERI+5], 1/20., spinor_field[DUM_DERI+5], VOLUME/2); 


  /* \delta Q sandwitched by X_e'^\dagger and X_o (EO)*/
  deriv_Sb(EO, DUM_DERI+4, DUM_DERI);
  deriv_Sb(EO, DUM_DERI+5, DUM_DERI+1);

  /* \delta Q sandwitched by X_o^\dagger and X_e  (OE)*/
  deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
  deriv_Sb(OE, DUM_DERI+1, DUM_DERI+3); 

}

void fermion_momenta_nond(double step) {
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
      deriv=&df0[i][mu];
      /* Factor 2 around? */
      tmp = 1.*step;
      _minus_const_times_mom(*xm,tmp,*deriv); 
    }
   }
}

