/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "linsolve.h"
#include "linalg_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "reweighting_factor.h"

double reweighting_factor(const int N, const double mu1, const double mu2) {
  int i, n_iter;
  double sq_norm, corr, sum=0.;

  /* Use spinor_field 2,3,5                         */
  /* in order not to conflict with anything else... */

  for(i = 0; i < N; i++) {
    random_spinor_field(2);
    g_mu = mu1;
    n_iter = solve_cg(3, 2, 0., 1.e-15, 1);
    
    g_mu = mu2;
    Qtm_pm_psi(spinor_field[5] , spinor_field[3]);
    
    sq_norm = square_norm(spinor_field[2], VOLUME/2);
    corr = scalar_prod_r(spinor_field[2], spinor_field[5], VOLUME/2);
    
    sq_norm -= corr;
    sum += sq_norm;
  }

  return(sum);
}

