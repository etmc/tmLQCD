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
#include"invert_eo.h"

int invert_eo(const int Even_new, const int Odd_new, const int Even, const int Odd) {

  int iter = 0;

  assign_mul_one_pm_imu_inv(Even_new, Even, +1.);
  
  Hopping_Matrix(OE, DUM_DERI, Even_new);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(DUM_DERI, 1., Odd, VOLUME/2);

  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  /* The solver inverts gamma_5 D ...          */
  gamma5(DUM_DERI, DUM_DERI); 
/*   iter = bicgstabell(Odd_new, DUM_DERI, 2000, 1.e-30, 2, 0.);   */
  iter = bicg(Odd_new, DUM_DERI, 0., 1.e-15);   

  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, DUM_DERI, Odd_new);
  mul_one_pm_imu_inv(DUM_DERI, +1.);
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new, 1., DUM_DERI, VOLUME/2); 

  return(iter);
}

void M_full(const int Even_new, const int Odd_new, const int Even, const int Odd) {
  /* Even sites */
  Hopping_Matrix(EO, DUM_DERI, Odd);
  assign_mul_one_pm_imu(Even_new, Even, 1.); 
  assign_add_mul_r(Even_new, -1., DUM_DERI, VOLUME/2);

  /* Odd sites */
  Hopping_Matrix(OE, DUM_DERI, Even);
  assign_mul_one_pm_imu(Odd_new, Odd, 1.); 
  assign_add_mul_r(Odd_new, -1., DUM_DERI, VOLUME/2);
}
