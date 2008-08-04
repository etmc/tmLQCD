/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "linsolve.h"
#include "linalg_eo.h"
#include "start.h"
#include "boundary.h"
#include "tm_operators.h"
#include "reweight_kappac.h"

double reweight_kappac(const int N) {
  int i, j, n_iter;
  double sq_norm, sq_norm_orig, corr, sum=0., sq_sum = 0., temp1;
  double mu_orig, kappa_orig, kappa_new;

  /* This routine is intended to do calculate reweighting factors for the 
     reweighting in kappa.
     For this purpose we first invert the TM-Dirac-operator at the kappa 
     as given in the input file and then multiply with a set of new TM-Dirac-
     operators with different kappa-values.
     In this way we'd like to check how far away from the original kappa the 
     reweighting remains reliable. */
  
  /* Use spinor_field 2,3,5                         */
  /* in order not to conflict with anything else... */

  /* First generate a random source: */
  random_spinor_field(g_spinor_field[2],VOLUME/2, 1);
  sq_norm_orig = square_norm(g_spinor_field[2], VOLUME/2);


  /* The invert with the original kappa and store the solution: 
     ----------------------------------------------------------*/
  /* store the original kappa, mu-value */
  kappa_orig = g_kappa;
  mu_orig = g_mu;
      
  /* This is the initial solution: */
  zero_spinor_field(g_spinor_field[3],VOLUME/2);

  /* and now solve:*/
  n_iter = solve_cg(g_spinor_field[3], g_spinor_field[2], 1.e-15, 1);




  /* Now run through the various new values of kappa and apply the 
     corresponding operator: */
  for(j = -10; j < 11; j++){
      kappa_new = kappa_orig + j*0.00001;
      /* the value of kappa is worked into the boundary conditions,
       so we need to update them: */
      boundary();
      
      /* Note that mu is in fact 2 kappa mu, so we have 
         2 mu kappa_new = 2 mu kappa_old * kappa_new/kappa_old */
      g_mu = mu_orig*kappa_new/kappa_orig;
      

      /* Multiply with the new operator: */
      Qtm_pm_psi(g_spinor_field[5] , g_spinor_field[3]);


      corr = scalar_prod_r(g_spinor_field[2], g_spinor_field[5], VOLUME/2);
      
      sq_norm = sq_norm_orig - corr;
      temp1 = sq_norm;
      sum += temp1;
      sq_sum += temp1*temp1;
      /*      
	      printf("rew: n_iter = %d, sq_norm = %e, corr = %e\n", n_iter, sq_norm, corr);
      */    
      printf("rew: kappa = %e, factor = %e, rew = %e\n", kappa_new, 2*sq_norm/VOLUME,exp(2*sq_norm/VOLUME));

  } /* end of loop over various new kappa-values. */

  return(sum);
}

