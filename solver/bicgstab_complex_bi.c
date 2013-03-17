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

  double err, squarenorm;
  _Complex double rho0, rho1, omega, alpha, beta, nom, denom;
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
  diff((spinor*)p, (spinor*)Q, (spinor*)r, 2*N);
  assign((spinor*)r, (spinor*)p, 2*N);
  assign((spinor*)hatr, (spinor*)p, 2*N);
  rho0 = scalar_prod((spinor*)hatr, (spinor*)r, 2*N, 1);
  squarenorm = square_norm((spinor*)Q, 2*N, 1);

  for(i = 0; i < max_iter; i++){
    err = square_norm((spinor*)r, 2*N, 1);
    if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
      printf("%d %e\n", i, err);
      fflush(stdout);
    }
  
    if((((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))) && i>0) {
      finalize_bisolver(bisolver_field, nr_sf);
      return(i);
    }
    f(v, p);
    denom = scalar_prod((spinor*)hatr, (spinor*)v, 2*N, 1);
    alpha = rho0 / denom;
    assign((spinor*)s, (spinor*)r, 2*N);
    assign_diff_mul((spinor*)s, (spinor*)v, alpha, 2*N);
    f(t, s);
    omega = scalar_prod((spinor*)t, (spinor*)s, 2*N, 1);
    omega /= square_norm((spinor*)t, 2*N, 1);
    assign_add_mul_add_mul((spinor*)P, (spinor*)p, (spinor*)s, alpha, omega, 2*N);
    assign((spinor*)r, (spinor*)s, 2*N);
    assign_diff_mul((spinor*)r, (spinor*)t, omega, 2*N);
    rho1 = scalar_prod((spinor*)hatr, (spinor*)r, 2*N, 1);
    nom = alpha * rho1;
    denom = omega * rho0;
    beta = nom / denom;
    omega = -omega;
    assign_mul_bra_add_mul_ket_add((spinor*)p, (spinor*)v, (spinor*)r, omega, beta, 2*N);
    rho0 = rho1;
  }
  finalize_bisolver(bisolver_field, nr_sf);
  return -1;
}
