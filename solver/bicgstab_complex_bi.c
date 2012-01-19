/***********************************************************************
 *
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
 *
 * The externally accessible functions are
 *
 *   int bicgstab(bispinor * const, bispinor * const, const int, double, matrix_mult_bi)
 *     BiCGstab solver
 * 
 *
 *
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 * 
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver_field.h"
#include "bicgstab_complex_bi.h"

/* P inout (guess for the solving bispinor)
   Q input
*/
int bicgstab_complex_bi(bispinor * const P, bispinor * const Q, const int max_iter, double eps_sq, const int rel_prec, const int N, matrix_mult_bi f){

  double err, d1, squarenorm;
  complex rho0, rho1, omega, alpha, beta, nom, denom;
  int i;
  bispinor * r, * p, * v, *hatr, * s, * t;
  bispinor ** bisolver_field = NULL;
  const int nr_sf = 6;

  if(N == VOLUME) {
    init_bisolver_field(&bisolver_field, VOLUMEPLUSRAND, nr_sf);
  }
  else {
    init_bisolver_field(&bisolver_field, VOLUMEPLUSRAND/2, nr_sf);
  }

  hatr = bisolver_field[0];
  r = bisolver_field[1];
  v = bisolver_field[2];
  p = bisolver_field[3];
  s = bisolver_field[4];
  t = bisolver_field[5];

  f(r, P);
  diff_bi(p, Q, r, N);
  assign_bi(r, p, N);
  assign_bi(hatr, p, N);
  rho0 = scalar_prod_bi(hatr, r, N);
  squarenorm = square_norm_bi(Q, N);

  for(i = 0; i < max_iter; i++){
    err = square_norm_bi(r, N);
    if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
      printf("%d %e\n", i, err);
      fflush(stdout);
    }
  
    if((((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) && i>0) {
      finalize_bisolver(bisolver_field, nr_sf);
      return(i);
    }
    f(v, p);
    denom = scalar_prod_bi(hatr, v, N);
    _div_complex(alpha, rho0, denom);
    assign_bi(s, r, N);
    assign_diff_mul_bi(s, v, alpha, N);
    f(t, s);
    omega = scalar_prod_bi(t,s, N);
    d1 = square_norm_bi(t, N);
    omega.re/=d1; omega.im/=d1;
    assign_add_mul_add_mul_bi(P, p, s, alpha, omega, N);
    assign_bi(r, s, N);
    assign_diff_mul_bi(r, t, omega, N);
    rho1 = scalar_prod_bi(hatr, r, N);
    _mult_assign_complex(nom, alpha, rho1);
    _mult_assign_complex(denom, omega, rho0);
    _div_complex(beta, nom, denom);
    omega.re=-omega.re; omega.im=-omega.im;
    assign_mul_bra_add_mul_ket_add_bi(p, v, r, omega, beta, N);
    rho0.re = rho1.re; rho0.im = rho1.im;
  }
  finalize_bisolver(bisolver_field, nr_sf);
  return -1;
}
