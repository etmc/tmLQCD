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

#include<stdlib.h>
#include"global.h"
#include"linalg_eo.h"
#include"tm_operators.h"
#include"Hopping_Matrix.h"
#include"linsolve.h"
#include"bicgstabell.h"
#include"clover_eo.h"
#include"solver/solver.h"
#include"invert_eo.h"

int invert_eo(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd,
	      const double precision, const int max_iter,
	      const int solver_flag) {

  int iter = 0;

  assign_mul_one_pm_imu_inv(Even_new, Even, +1.);
  
  Hopping_Matrix(OE, spinor_field[DUM_DERI], Even_new);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(spinor_field[DUM_DERI], +1., Odd, VOLUME/2);

  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  /* The solver inverts gamma_5 D ...          */
  gamma5(DUM_DERI, DUM_DERI); 
  /*   iter = bicgstabell(Odd_new, DUM_DERI, 2000, 1.e-30, 2, 0.);   */
  /*   iter = bicg(Odd_new, spinor_field[DUM_DERI], 0., 1.e-15);    */

  if(solver_flag == BICGSTAB) {
    if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
    iter = bicgstab_complex(Odd_new, spinor_field[DUM_DERI], max_iter, precision, &Qtm_plus_psi);
  }
  else if(solver_flag == GMRES) {
    if(g_proc_id == 0) {printf("# Using GMRES!\n"); fflush(stdout);}
    iter = gmres(Odd_new, spinor_field[DUM_DERI], 10, max_iter/10, precision, &Qtm_plus_psi);
  }
  else if(solver_flag == CG) {
    if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
    iter = cg_her(Odd_new, spinor_field[DUM_DERI], max_iter, precision, &Qtm_pm_psi, 0, 0.);
    Qtm_minus_psi(Odd_new, Odd_new);
  }
  if(solver_flag == MR) {
    if(g_proc_id == 0) {printf("# Using MR!\n"); fflush(stdout);}
    iter = mr(Odd_new, spinor_field[DUM_DERI], max_iter, precision, &Qtm_plus_psi);
  }
  else if(solver_flag == CGS) {
    if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}
    iter = cgs_real(Odd_new, spinor_field[DUM_DERI], max_iter, precision, &Qtm_plus_psi);
  }

  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, spinor_field[DUM_DERI], Odd_new);
  mul_one_pm_imu_inv(spinor_field[DUM_DERI], +1.);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new, spinor_field[DUM_DERI], +1., VOLUME/2); 
  
  return(iter);
}

void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, spinor_field[DUM_DERI], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1.); 
  assign_add_mul_r(Even_new, spinor_field[DUM_DERI], -1., VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, spinor_field[DUM_DERI], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1.); 
  assign_add_mul_r(Odd_new, spinor_field[DUM_DERI], -1., VOLUME/2);
}
