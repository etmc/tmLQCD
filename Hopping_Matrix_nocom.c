/* $Id$ */

/******************************************
 * Hopping_Matrix is the conventional Wilson 
 * hopping matrix
 *
 *
 * But the communication is left out by a 
 * dirty trick...
 *
 * \kappa\sum_{\pm\mu}(r+\gamma_\mu)U_{x,\mu}
 *
 * for ieo = 0 this is M_{eo}, for ieo = 1
 * it is M_{oe}
 *
 * l is the number of the output field
 * k is the number of the input field
 *
 ******************************************/

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "xchange.h"
#include "su3.h"
#include "sse.h"
#include "boundary.h"
#include "Hopping_Matrix.h"

#define xchange_field(a) 
#define Hopping_Matrix Hopping_Matrix_nocom

#include "Hopping_Matrix.c"
