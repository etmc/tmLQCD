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
#include "clover_eo.h" 
#include "start.h"
#include "sw.h"
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
  double qo;

  for(i=0;i<(VOLUME+RAND);i++){ 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }

  /* If CG is used anyhow */
  /* don't know, why g_5 multiplication 
     gamma5(spinor_field[DUM_DERI+1], spinor_field[nond_psf_strange], VOLUME/2);
     gamma5(spinor_field[DUM_DERI+2], spinor_field[nond_psf_charm], VOLUME/2);
  */
  /* Invert Q_{+} Q_{-} */
  /* X_o -> DUM_DERI+1 */

  /* This replaces the solver */
  QdaggerQ_power(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2], 
		 dop_cheby_coef, dop_n_cheby, 
		 spinor_field[STRANGE], spinor_field[CHARM]);
  
  /* Y_o -> DUM_DERI  */
  /* check this replacement */
/*   QNon_degenerate(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI+4],spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2]); */
  M_full(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI+4],spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2]);

  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X */
  /* this is checked with mubar and epsbar=0 apart from  */
  /* check whether mul_one_minus_imubar of mul_one_minus_imubar_inv */
  QNon_degenerate_eo(spinor_field[DUM_DERI+5], spinor_field[DUM_DERI+6], 
		     spinor_field[DUM_DERI+3], spinor_field[DUM_DERI+4]); 
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, DUM_DERI+1, DUM_DERI+5);
  deriv_Sb(EO, DUM_DERI+6, DUM_DERI+2);

  /* to get the even sites of Y */
  /* this is checked with mubar and epsbar=0 apart from  */
  /* check whether mul_one_minus_imubar of mul_one_minus_imubar_inv */
  QNon_degenerate_eo(spinor_field[DUM_DERI+5], spinor_field[DUM_DERI+6], 
		     spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2]); 
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(OE, DUM_DERI+3, DUM_DERI+5);
  deriv_Sb(EO, DUM_DERI+6, DUM_DERI+4);

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
      /* This 2* is coming from what?             */
      /* From a missing factor 2 in trace_lambda? */
      tmp = 2.*step;
      _minus_const_times_mom(*xm,tmp,*deriv); 
    }
   }
}

