/* $Id$ */

/************************************************
 *
 * Typedefinition of the pointer to the function
 * which contains the matrix multiplication.
 *
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 ************************************************/

#ifndef _MATRIX_MULT_TYPEDEF_BI_H
#define _MATRIX_MULT_TYPEDEF_BI_H

typedef void (*matrix_mult_bi) (bispinor * const, bispinor * const);

/*   SO FAR, THE CLOVER TERM IS NOT REALLY NEEDED FOR THE ND-CASE */
/*
typedef void (*matrix_mult_clover) (spinor * const, spinor * const, const double);
*/

#endif
