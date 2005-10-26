/**************************************************************************
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver/matrix_mult_typedef.h"
#include "cgs_real.h"

#ifdef _SOLVER_OUTPUT
#define _SO(x) x
#else
#define _SO(x)
#endif

/* P inout (guess for the solving spinor)
   Q input
*/

int cgs_real(spinor * const P, spinor * const Q, const int max_iter, 
	     double eps_sq, const int rel_prec, const int N, matrix_mult f) {
  static complex alpha, beta,rjr0,nom,denom,one;
  static double res_sq, squarenorm;
  int i;

/*   init_solver_field(6); */
  one.re=1.;one.im=0.;
  /* Initialisierung der sf-Felder */  
  f(spinor_field[DUM_SOLVER],P);
  diff(spinor_field[DUM_SOLVER],Q,spinor_field[DUM_SOLVER], N); /* residual in sf0 */
  assign(spinor_field[DUM_SOLVER+1],spinor_field[DUM_SOLVER], N); 
  assign(spinor_field[DUM_SOLVER+2],spinor_field[DUM_SOLVER], N);
  assign(spinor_field[DUM_SOLVER+5],spinor_field[DUM_SOLVER], N); /* ri=pi=ui=r0 */
  squarenorm = square_norm(Q, N);

  /* loop! */
  for(i=0;i<=max_iter;i++) {
    res_sq=square_norm(spinor_field[DUM_SOLVER], N);
    if(g_proc_id == g_stdio_proc) {
      printf("%d\t%g\n",i,res_sq); 
      fflush( stdout );
    }
    rjr0.re=scalar_prod_r(spinor_field[DUM_SOLVER],spinor_field[DUM_SOLVER+5], N);
    /*     square_and_prod(&res_sq,&rjr0,spinor_field[DUM_SOLVER],spinor_field[DUM_SOLVER+5]); */
    if(((res_sq<eps_sq) && (rel_prec == 0)) || ((res_sq<eps_sq*squarenorm) && (rel_prec == 1))) {
      return i;
    }
    f(spinor_field[DUM_SOLVER+3],spinor_field[DUM_SOLVER+1]);	/* calc v */
    /* calc alpha */
    denom.re=scalar_prod_r(spinor_field[DUM_SOLVER+3],spinor_field[DUM_SOLVER+5], N);
    /* _div_complex(alpha,rjr0,denom);*/
    alpha.re=rjr0.re/denom.re;
    /* calc q */
    mul_diff_mul(spinor_field[DUM_SOLVER+4],spinor_field[DUM_SOLVER+2],spinor_field[DUM_SOLVER+3],one,alpha, N); 
    /* calc P and residual */
    /* calc alpha(u+q) into sf2 */
    assign_add_mul_r(spinor_field[DUM_SOLVER+2],spinor_field[DUM_SOLVER+4],1., N);
    assign_add_mul(P,spinor_field[DUM_SOLVER+2],alpha, N); /* new P */
    /* calc new residual */
    f(spinor_field[DUM_SOLVER+3],spinor_field[DUM_SOLVER+2]);
    assign_diff_mul(spinor_field[DUM_SOLVER],spinor_field[DUM_SOLVER+3],alpha, N);
    /* calc beta */
    nom.re=scalar_prod_r(spinor_field[DUM_SOLVER],spinor_field[DUM_SOLVER+5], N);
    /* _div_complex(beta,nom,rjr0); */
    beta.re=nom.re/rjr0.re;
    /* calc new u */
    mul_add_mul(spinor_field[DUM_SOLVER+2],spinor_field[DUM_SOLVER],spinor_field[DUM_SOLVER+4],one,beta, N);
    /* calc new p */
    /* _mult_assign_complex(nom,beta,beta); */
    nom.re=beta.re*beta.re;
    mul(spinor_field[DUM_SOLVER+1],nom,spinor_field[DUM_SOLVER+1], N);	
    assign_add_mul_add_mul(spinor_field[DUM_SOLVER+1],spinor_field[DUM_SOLVER+4],spinor_field[DUM_SOLVER+2],beta,one, N);
  }
  return -1;
}

/*
mul_diff_mul
mul

*/





