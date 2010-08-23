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
/**************************************************************************
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
#include "linalg/mul_diff_mul_r.h"
#include "linalg/assign_add_mul_add_mul_r.h"
#include "solver/matrix_mult_typedef.h"
#include "cgs_real.h"


/* P inout (guess for the solving spinor)
   Q input
*/

int cgs_real(spinor * const P, spinor * const Q, const int max_iter, 
	     double eps_sq, const int rel_prec, const int N, matrix_mult f) {
  static double alpha, beta,rjr0,nom,denom,one;
  static double res_sq, squarenorm;
  int i;

/*   init_solver_field(6); */
  one=1.;
  /* Initialisierung der sf-Felder */  
  f(g_spinor_field[DUM_SOLVER],P);
  diff(g_spinor_field[DUM_SOLVER],Q,g_spinor_field[DUM_SOLVER], N); /* residual in sf0 */
  assign(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER], N); 
  assign(g_spinor_field[DUM_SOLVER+2],g_spinor_field[DUM_SOLVER], N);
  assign(g_spinor_field[DUM_SOLVER+5],g_spinor_field[DUM_SOLVER], N); /* ri=pi=ui=r0 */
  squarenorm = square_norm(Q, N, 1);

  /* loop! */
  for(i=0;i<=max_iter;i++) {
    res_sq=square_norm(g_spinor_field[DUM_SOLVER], N, 1);
    if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
      printf("%d\t%g\n",i,res_sq); 
      fflush( stdout );
    }
    rjr0 = scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+5], N, 1);
    /*     square_and_prod(&res_sq,&rjr0,g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER+5]); */
    if(((res_sq<eps_sq) && (rel_prec == 0)) || ((res_sq<eps_sq*squarenorm) && (rel_prec == 1))) {
      return i;
    }
    f(g_spinor_field[DUM_SOLVER+3],g_spinor_field[DUM_SOLVER+1]);	/* calc v */
    /* calc alpha */
    denom=scalar_prod_r(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+5], N, 1);
    /* _div_complex(alpha,rjr0,denom);*/
    alpha=rjr0/denom;
    /* calc q */
    mul_diff_mul_r(g_spinor_field[DUM_SOLVER+4], g_spinor_field[DUM_SOLVER+2], 
		   g_spinor_field[DUM_SOLVER+3],one,alpha, N); 
    /* calc P and residual */
    /* calc alpha(u+q) into sf2 */
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+4],1., N);
    assign_add_mul_r(P,g_spinor_field[DUM_SOLVER+2], alpha, N); /* new P */
    /* calc new residual */
    f(g_spinor_field[DUM_SOLVER+3],g_spinor_field[DUM_SOLVER+2]);
    assign_add_mul_r(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+3], -alpha, N);
    /* calc beta */
    nom=scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+5], N, 1);
    /* _div_complex(beta,nom,rjr0); */
    beta = nom/rjr0;
    /* calc new u */
    mul_diff_mul_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], 
		   g_spinor_field[DUM_SOLVER+4], one, -beta, N);
    /* calc new p */
    /* _mult_assign_complex(nom,beta,beta); */
    nom=beta*beta;
    mul_r(g_spinor_field[DUM_SOLVER+1], nom, g_spinor_field[DUM_SOLVER+1], N);	
    assign_add_mul_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+4],
			     g_spinor_field[DUM_SOLVER+2], beta, one, N);
  }
  return -1;
}

/*
mul_diff_mul
mul

*/





