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

/*******************************************************************************
 * Generalized minimal residual (GMRES) with a maximal number of restarts.    
 * Solves Q=AP for _Complex double regular matrices A.
 * For details see: Andreas Meister, Numerik linearer Gleichungssysteme        
 *   or the original citation:                                                 
 * Y. Saad, M.H.Schultz in GMRES: A generalized minimal residual algorithm    
 *                         for solving nonsymmetric linear systems.            
 * 			SIAM J. Sci. Stat. Comput., 7: 856-869, 1986           
 *           
 * int gmres(spinor * const P,spinor * const Q, 
 *	   const int m, const int max_restarts,
 *	   const double eps_sq, matrix_mult f)
 *                                                                 
 * Returns the number of iterations needed or -1 if maximal number of restarts  
 * has been reached.                                                           
 *
 * Inout:                                                                      
 *  spinor * P       : guess for the solving spinor                                             
 * Input:                                                                      
 *  spinor * Q       : source spinor
 *  int m            : Maximal dimension of Krylov subspace                                     
 *  int max_restarts : maximal number of restarts                                   
 *  double eps       : stopping criterium                                                     
 *  matrix_mult f    : pointer to a function containing the matrix mult
 *                     for type matrix_mult see matrix_mult_typedef.h
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 ********************************************************************************/

#ifndef _GMRES_DR_H
#define _GMRES_DR_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int gmres_dr(spinor * const P,spinor * const Q, 
	     const int m, const int nr_ev, const int max_restarts,
	     const double eps, const int rel_prec, 
	     const int N, matrix_mult f);

#endif
