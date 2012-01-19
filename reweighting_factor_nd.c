/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#include "linsolve.h"
#include "linalg_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "Ptilde_nd.h"
#include "phmc.h"
#include "reweighting_factor_nd.h"

double reweighting_factor_nd(const int N) {
  int i, n_iter;
  double sq_norm, corr, sum=0., sq_sum = 0., temp1;
  double mu1, mu2;

  complex temp2;

  mu1 = g_mu1;
  mu2 = g_mu1;
  
  /* Use spinor_field 2,3,5                         */
  /* in order not to conflict with anything else... */

  for(i = 0; i < N; i++) {
    random_spinor_field(g_chi_up_spinor_field[2],VOLUME/2, 1);
    random_spinor_field(g_chi_dn_spinor_field[2],VOLUME/2, 1);
    zero_spinor_field(g_chi_up_spinor_field[3],VOLUME/2);
    zero_spinor_field(g_chi_dn_spinor_field[3],VOLUME/2);

    temp1 = phmc_ptilde_cheby_coef[0];
    phmc_ptilde_cheby_coef[0] = temp1 - 1;

    Poly_tilde_ND(g_chi_up_spinor_field[3], g_chi_dn_spinor_field[3], phmc_ptilde_cheby_coef, phmc_ptilde_n_cheby, g_chi_up_spinor_field[2], g_chi_dn_spinor_field[2]);

    phmc_ptilde_cheby_coef[0] = temp1;

    temp2 = scalar_prod(g_chi_up_spinor_field[2], g_chi_up_spinor_field[3], VOLUME/2, 1);
    if(temp2.im > 1.0e-8) {
      printf("!!! WARNING  Immaginary part of CORR-UP  LARGER than 10^-8 !!! \n");
      printf(" CORR-UP:  Re=%12.10e  Im=%12.10e \n", temp2.re, temp2.im);
    }
    corr = temp2.re;
    printf(" CORR-UP:  Re=%12.10e \n", corr);
    temp2 = scalar_prod(g_chi_dn_spinor_field[2], g_chi_dn_spinor_field[3], VOLUME/2, 1);
    if(temp2.im > 1.0e-8) {
      printf("!!! WARNING  Immaginary part of CORR_DN  LARGER than 10^-8 !!! \n");
      printf(" CORR-DN:  Re=%12.10e  Im=%12.10e \n", temp2.re, temp2.im);
    }
    corr += temp2.re;
    printf(" CORR-DN:  Re=%12.10e \n", temp2.im);

    temp1 = -corr;
    sum += temp1;
    sq_sum += temp1*temp1;
    printf("rew: n_iter = %d, sq_norm = %e, corr = %e\n", n_iter, sq_norm, corr);

    /*    
    random_spinor_field(g_spinor_field[2],VOLUME/2, 1);
    g_mu = mu2;
    zero_spinor_field(g_spinor_field[3],VOLUME/2);
    n_iter = solve_cg(3, 2, 0., 1.e-15, 1);

    g_mu = mu1;
    Qtm_pm_psi(g_spinor_field[5] , g_spinor_field[3]);

    sq_norm = square_norm(g_spinor_field[2], VOLUME/2, 1);
    corr = scalar_prod_r(g_spinor_field[2], g_spinor_field[5], VOLUME/2, 1);
    
    sq_norm -= corr;
    temp1 = sq_norm;
    sum += temp1;
    sq_sum += temp1*temp1;
    printf("rew: n_iter = %d, sq_norm = %e, corr = %e\n", n_iter, sq_norm, corr);
    */

  }
  sum/=(double)N;
  sq_sum/=(double)N;
  printf("rew: factor = %e, err = %e\n", sum, sqrt(sum*sum-sq_sum)/((double)N-1));
  return(sum);
}

