/* $Id$ */

/****************************************************************
 *
 * invert_eo makes an inversion with EO precoditioned
 * tm Operator
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_eo returns the number of iterations neede or -1 if the 
 * solver did not converge.
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ****************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include"global.h"
#include"linalg_eo.h"
#include"tm_operators.h"
#include"Hopping_Matrix.h"
#include"D_psi.h"
#include"linsolve.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange.h"
#include"Nondegenerate_Matrix.h"
#include"invert_doublet_eo.h"

int invert_doublet_eo(spinor * const Even_new_s, spinor * const Even_new_c, 
		      spinor * const Odd_new_c, spinor * const Odd_new_s, 
		      spinor * const Even_s, spinor * const Even_c,
		      spinor * const Odd_s, spinor * const Odd_c,
		      const double precision, const int max_iter,
		      const int solver_flag, const int rel_prec) {

  int iter = 0;

  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using Even/Odd preconditioning!\n"); fflush(stdout);}
  M_ee_inv_ND(Even_new_s, Even_new_c, 
	      Even_s, Even_c);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new_s);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI+1], Even_new_c);

  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd_s, VOLUME/2);
  assign_mul_add_r(g_spinor_field[DUM_DERI+1], +1., Odd_c, VOLUME/2);
/*   assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2); */
  
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  
  /* Here we invert the hermitean operator squared */
  if(g_proc_id == 0) {
    printf("# Using CG for flavour doublet!\n"); 
    fflush(stdout);
  }
  iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		   max_iter, precision, rel_prec, 
		   VOLUME/2, &Q_Qdagger_ND, 0, 1000);
  QdaggerNon_degenerate(Odd_new_s, Odd_new_c,
			Odd_new_s, Odd_new_c);
/*   Qtm_minus_psi(Odd_new, Odd_new); */
  
  
  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  M_ee_inv_ND(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
	      g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
/*   assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2); */

  return(iter);
}

