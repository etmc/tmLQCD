/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "linsolve.h"
#include "deriv_Sb.h"
#include "clover_eo.h"
#include "tm_operators.h"
#include "hybrid_update.h"
#include "Hopping_Matrix.h"
#include "derivative_psf.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;

void derivative_psf(const int nr) {

  int i,mu;
  double sum=0.;
  su3adj * deriv;

  for(i=0;i<(VOLUME+RAND);i++){ 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }

  if(nr == 0) {
    g_mu = g_mu1;
    if(ITER_MAX_BCG == 0){
      /* If CG is used anyhow */
      gamma5(DUM_DERI+1, first_psf);
      /* Invert Q_{+} Q_{-} */
      /* X_o -> DUM_DERI+1 */
      count00 += solve_cg(DUM_DERI+1, first_psf, 0., EPS_SQ1);
      /* Y_o -> DUM_DERI  */
      Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
    }
    else{
      /*contributions from field 0 -> first_psf*/
      gamma5(DUM_DERI, first_psf);
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      count00 += bicg(DUM_DERI, first_psf, 0., EPS_SQ1);
      gamma5(DUM_DERI+1,DUM_DERI);
      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      count01 += bicg(DUM_DERI+1,DUM_DERI,0.,EPS_SQ1);
      g_mu = -g_mu;
    }
  }
  else if(nr == 1) {
    /* First term coming from the second field */
    /* Multiply with W_+ */
    g_mu = g_mu1;	
    Qtm_plus_psi(spinor_field[DUM_DERI+2], spinor_field[second_psf]);
    g_mu = g_mu2;
    if(ITER_MAX_BCG == 0){
      /* If CG is used anyhow */
      gamma5(DUM_DERI+1, DUM_DERI+2);
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      count10 += solve_cg(DUM_DERI+1, DUM_DERI+2, 0., EPS_SQ1);
      /* Y_W -> DUM_DERI  */
      Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
    }
    else{
      gamma5(DUM_DERI, DUM_DERI+2);
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      count10 += bicg(DUM_DERI, DUM_DERI+2, 0., EPS_SQ1);
      gamma5(DUM_DERI+1,DUM_DERI);
      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      count11 += bicg(DUM_DERI+1,DUM_DERI,0.,EPS_SQ1);
      g_mu = -g_mu;   
    }

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+2], spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
    g_mu = g_mu1;
    
    
    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(spinor_field[DUM_DERI], -1., spinor_field[second_psf], VOLUME/2);
    g_mu = g_mu1;
  }

  if(nr == 2) {
    /* First term coming from the third field */
    /* Multiply with W_+ */
    g_mu = g_mu2;	
    Qtm_plus_psi(spinor_field[DUM_DERI+2], spinor_field[third_psf]);
    g_mu = g_mu3;
    if(ITER_MAX_BCG == 0){
      /* If CG is used anyhow */
      gamma5(DUM_DERI+1, DUM_DERI+2);
      /* Invert Q_{+} Q_{-} */
      /* X_W -> DUM_DERI+1 */
      count20 += solve_cg(DUM_DERI+1, DUM_DERI+2, 0., EPS_SQ1);
      /* Y_W -> DUM_DERI  */
      Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
    }
    else{
      gamma5(DUM_DERI, DUM_DERI+2);
      /* Invert first Q_+ */
      /* Y_o -> DUM_DERI  */
      count20 += bicg(DUM_DERI, DUM_DERI+2, 0., EPS_SQ1);
      gamma5(DUM_DERI+1,DUM_DERI);
      /* Now Q_- */
      /* X_o -> DUM_DERI+1 */
      g_mu = -g_mu;
      count21 += bicg(DUM_DERI+1,DUM_DERI,0.,EPS_SQ1);
      g_mu = -g_mu;   
    }

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+2], spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
    g_mu = g_mu1;

    /* Second term coming from the third field */
    /* The sign is opposite!! */
    mul_r( spinor_field[DUM_DERI], -1., spinor_field[third_psf], VOLUME/2);
    g_mu = g_mu2;
  }

  /* apply Hopping Matrix M_{eo} */
  /* to get the even sites of X */
  H_eo_tm_inv_psi(spinor_field[DUM_DERI+2], spinor_field[DUM_DERI+1], EO, -1.);
  /* \delta Q sandwitched by Y_o^\dagger and X_e */
  deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
  
  /* to get the even sites of Y */
  H_eo_tm_inv_psi(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI], EO, +1);
  /* \delta Q sandwitched by Y_e^\dagger and X_o */
  deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
  g_mu = g_mu1;
  
  for(i = 0; i < VOLUME; i++){
    for(mu=0;mu<4;mu++){
      deriv=&df0[i][mu];
      sum+= _su3adj_square_norm(*deriv);
    }
  }
  if(g_proc_id == 0) {
    printf("nach %d fermiongesamt %e\n", nr, sum/((double)VOLUME));
    fflush(stdout);
    sum=0.;
  }
}
