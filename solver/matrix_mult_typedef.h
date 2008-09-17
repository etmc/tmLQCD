/* $Id$ */

/************************************************
 *
 * Typedefinition of the pointer to the function
 * which contains the matrix multiplication.
 *
 ************************************************/

#ifndef _MATRIX_MULT_TYPEDEF_H
#define _MATRIX_MULT_TYPEDEF_H

typedef void (*matrix_mult) (spinor * const, spinor * const);
typedef void (*matrix_mult_clover) (spinor * const, spinor * const, const double);
typedef void (*c_matrix_mult) (complex * const, complex * const);

#endif
