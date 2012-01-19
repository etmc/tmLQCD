/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

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
       const int parallel, matrix_mult f);

int mrblk(spinor * const P, spinor * const Q,
	  const int max_iter, const double eps_sq,
	  const int rel_prec, const int N, 
	  matrix_mult_blk f, const int blk);

#endif
