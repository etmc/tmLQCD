/**************************************************************************
 *
 * $Id$
 *
 * The externally accessible functions are
 *
 *   int bicgstab(spinor * const, spinor * const, const int, double, matrix_mult)
 *     BiCGstab solver
 * 
 *
 *
 *
 * Author: Carsten Urbach 
 *         <urbach@ifh.de>
 * 
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "bicgstab_complex.h"

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif
 


/* P inout (guess for the solving spinor)
   Q input
*/
int bicgstab_complex(spinor * const P,spinor * const Q, const int max_iter, 
		     double eps_sq, matrix_mult f){
  double err, d1;
  complex rho0, rho1, omega, alpha, beta, nom, denom;
  int i, N=VOLUME/2;
  spinor * r, * p, * v, *hatr, * s, * t;

/*   init_solver_field(6); */
  hatr = spinor_field[DUM_SOLVER];
  r = spinor_field[DUM_SOLVER+1];
  v = spinor_field[DUM_SOLVER+2];
  p = spinor_field[DUM_SOLVER+3];
  s = spinor_field[DUM_SOLVER+4];
  t = spinor_field[DUM_SOLVER+5];

  f(r, P);
  diff(p, Q, r, N);
  assign(r, p, N);
  assign(hatr, p, N);
  rho0 = scalar_prod(hatr, r, N);

  for(i = 0; i < max_iter; i++){
    err = square_norm(r, N);
    _SO(if(g_proc_id == g_stdio_proc){printf("%d %e\n", i, err);} );
    _SO(if(g_proc_id == g_stdio_proc){fflush(stdout);} );
    if(err < eps_sq && i>0){
      return(i);
    }
    f(v, p);
    denom = scalar_prod(hatr, v, N);
    _div_complex(alpha, rho0, denom);
    assign(s, r, N);
    assign_diff_mul(s, v, alpha, N);
    f(t, s);
    omega = scalar_prod(t,s, N);
    d1 = square_norm(t, N);
    omega.re/=d1; omega.im/=d1;
    assign_add_mul_add_mul(P, p, s, alpha, omega, N);
    assign(r, s, N);
    assign_diff_mul(r, t, omega, N);
    rho1 = scalar_prod(hatr, r, N);
    _mult_assign_complex(nom, alpha, rho1);
    _mult_assign_complex(denom, omega, rho0);
    _div_complex(beta, nom, denom);
    omega.re=-omega.re; omega.im=-omega.im;
    assign_mul_bra_add_mul_ket_add(p, v, r, omega, beta, N);
    rho0.re = rho1.re; rho0.im = rho1.im;
  }
  return -1;
}
