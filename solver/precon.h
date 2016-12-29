/*********************************************************************** 
 * preconditioning related functions.
 * Note: there exist already polynomial preconditioning operators
 *       related to the chebychev polynomials in poly_precon.h files
 *       that was developed earlier by Carsten
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/

#ifndef _PRECON_H
#define _PRECON_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "linalg_eo.h"
/* #include "memalloc.h" */
#include "start.h"

void cheb_poly_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k, spinor *v1, spinor *v2);
/*
 Chebyshev polynomial of the oprator f in the interval [a,b] normalized to 1 at 0.
 R: output spinor
 S: input spinor
 f: matrix-vector multiplication for the operator
 a,b: limits of the interval with b>a
 k: degree of the polynomial
 v1 and v2 are spinors needed to be used in the recurrence relation
*/


void cheb_poly_roots(_Complex double *roots, const int k, const double a, const double b);
/*
roots of the shifted Chebyshev polynomial of degree k in the interval [a,b]
*/


#endif
