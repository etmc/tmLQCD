/* $Id$ */
/***********************************************
 * This function defines the boundary cond.
 * 
 *
 ***********************************************/

#include <stdlib.h>
#include "math.h"
#include "global.h"
#include "su3.h"
#include "boundary.h"

complex ka0, ka1, ka2, ka3;
const double PI_ = 3.14159265358979;
double X0;

void boundary(){
  double x0,x1,x2,x3;
  /* anti-periodic in time */
  x0 = X0 * PI_/((T)*g_nproc);
  x1 = X1 * PI_/(L);
  x2 = X2 * PI_/(L);
  x3 = X3 * PI_/(L);
  ka0.re = g_kappa * cos(x0); 
  ka0.im = g_kappa * sin(x0);
  ka1.re = g_kappa * cos(x1); 
  ka1.im = g_kappa * sin(x1);
  ka2.re = g_kappa * cos(x2); 
  ka2.im = g_kappa * sin(x2);
  ka3.re = g_kappa * cos(x3); 
  ka3.im = g_kappa * sin(x3);
}
