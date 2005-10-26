/* $Id$ */

/****************************************************
 * Minimal residual solver
 * int mr(spinor * const P, spinor * const Q,
 *	const int max_iter, const double eps_sq,
 *	matrix_mult f){ *
 *
 * returns the number of iterations needed to reach
 * the desired precision. return -1 if the maximal
 * number of iterations was reached.
 *
 * Inout:                                                                      
 *  spinor * P       : guess for the solving spinor                                             
 * Input:                                                                      
 *  spinor * Q       : source spinor
 *  int max_iter     : maximal number of iterations                                 
 *  double eps_sqr   : stopping criterium                                                     
 *  matrix_mult f    : pointer to a function containing 
 *                     the matrix mult for type 
 *                     matrix_mult see 
 *                     matrix_mult_typedef.h
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

#ifndef _MR_H
#define _MR_H

int mr(spinor * const P, spinor * const Q,
       const int max_iter, const double eps_sq,
       const int rel_prec, const int N, 
       matrix_mult f);

#endif
