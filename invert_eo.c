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
#include"linsolve.h"
#include"gamma.h"
#include"solver/solver.h"
#include"xchange.h"
#include"invert_eo.h"

int invert_eo(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd,
	      const double precision, const int max_iter,
	      const int solver_flag, const int rel_prec) {

  int iter = 0;

  assign_mul_one_pm_imu_inv(Even_new, Even, +1.);
  
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new); 
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd, VOLUME/2);

  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */

  if(solver_flag == BICGSTAB) {
    if(g_proc_id == 0) {printf("# Using BiCGstab!\n"); fflush(stdout);}
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
    iter = bicgstab_complex(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
  }
  else if(solver_flag == GMRES) {
    if(g_proc_id == 0) {printf("# Using GMRES!\n"); fflush(stdout);}
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
    iter = gmres(Odd_new, g_spinor_field[DUM_DERI], 10, max_iter/10, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
  }
  else if(solver_flag == FGMRES) {
    if(g_proc_id == 0) {printf("# Using GMRES!\n"); fflush(stdout);}
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.);
    iter = fgmres(Odd_new, g_spinor_field[DUM_DERI], 10, max_iter/10, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
  }
  else if(solver_flag == BICGSTABELL) {
    if(g_proc_id == 0) {printf("# Using BiCGstab2!\n"); fflush(stdout);}
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
    iter = bicgstabell(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, 3, VOLUME/2, &Mtm_plus_sym_psi);
    mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
  }
  else if(solver_flag == CG) {
    /* Here we invert the hermitean operator squared */
    gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
    if(g_proc_id == 0) {printf("# Using CG!\n"); fflush(stdout);}
    iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi, 0, 0.);
    Qtm_minus_psi(Odd_new, Odd_new);
  }
  else if(solver_flag == MR) {
    if(g_proc_id == 0) {printf("# Using MR!\n"); fflush(stdout);}
    iter = mr(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_psi);
  }
  else if(solver_flag == CGS) {
    if(g_proc_id == 0) {printf("# Using CGS!\n"); fflush(stdout);}
    mul_one_pm_imu_inv(g_spinor_field[DUM_DERI], +1.); 
    iter = cgs_real(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Mtm_plus_sym_psi);
    mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.); 
  }
  else {
    if(g_proc_id == 0) {printf("# Using CG as default solver!\n"); fflush(stdout);}
    gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
    iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi, 0, 0.);
    Qtm_minus_psi(Odd_new, Odd_new);
  }

  /* In case of failure, redo with CG */
  if(iter == -1 && solver_flag !=CG) {
    /* Here we invert the hermitean operator squared */
    gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);  
    if(g_proc_id == 0) {printf("# Redoing it with CG!\n"); fflush(stdout);}
    iter = cg_her(Odd_new, g_spinor_field[DUM_DERI], max_iter, precision, rel_prec, VOLUME/2, &Qtm_pm_psi, 0, 0.);
    Qtm_minus_psi(Odd_new, Odd_new);
  }

  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new);
  mul_one_pm_imu(g_spinor_field[DUM_DERI], +1.);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], +1., VOLUME/2); 
  
  return(iter);
}

void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd) {
  /* Even sites */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1.); 
  assign_add_mul_r(Even_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1.); 
  assign_add_mul_r(Odd_new, g_spinor_field[DUM_DERI], -1., VOLUME/2);
}
