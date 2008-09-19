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
#include "tm_operators.h"
#include "reweighting_factor.h"

double reweighting_factor(const int N) {
  int i, n_iter;
  double sq_norm, corr, sum=0., sq_sum = 0., temp1;
  double mu1, mu2;

  if(g_nr_of_psf == 2) {
    /* The physical mass */
    mu1 = g_mu2;
    /* The precond. mass */
    mu2 = g_mu1;
  }
  else if (g_nr_of_psf == 3) {
    /* The physical mass */
    mu1 = g_mu3;
    /* The precond. mass */
    mu2 = g_mu2;
  }

  /* Use spinor_field 2,3,5                         */
  /* in order not to conflict with anything else... */

  for(i = 0; i < N; i++) {
    random_spinor_field(g_spinor_field[2],VOLUME/2);
    g_mu = mu2;
    zero_spinor_field(g_spinor_field[3],VOLUME/2);
    n_iter = solve_cg(g_spinor_field[3], g_spinor_field[2], 0., 1.e-15, 1);
    
    g_mu = mu1;
    Qtm_pm_psi(g_spinor_field[5] , g_spinor_field[3]);
    
    sq_norm = square_norm(g_spinor_field[2], VOLUME/2, 1);
    corr = scalar_prod_r(g_spinor_field[2], g_spinor_field[5], VOLUME/2, 1);
    
    sq_norm -= corr;
    temp1 = sq_norm;
    sum += temp1;
    sq_sum += temp1*temp1;
    printf("rew: n_iter = %d, sq_norm = %e, corr = %e\n", n_iter, sq_norm, corr);
  }
  sum/=(double)N;
  sq_sum/=(double)N;
  printf("rew: factor = %e, err = %e\n", sum, sqrt(sum*sum-sq_sum)/((double)N-1));
  return(sum);
}

