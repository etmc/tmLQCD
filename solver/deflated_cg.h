/*****************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Deflating CG using eigenvectors computed using ARPACK
 *
 * Author: A.M. Abdel-Rehim (amabdelrehim@gmail.com)
 *         Novemebr 2014
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
 *
 ****************************************************************************/
/* A sample input is given in the sample-input folder */

#ifndef _DEFLATED_CG_H
#define _DEFLATED_CG_H

#include "su3.h"
/* #include "solver/matrix_mult_typedef.h" */
#include "solver/solver_params.h"
#include "solver/deflator.h"


int exactdeflated_cg(
  solver_params_t solver_params, /* (IN) parameters for solver */
  deflator_params_t deflator_params, /* (IN) parameters for deflator */
  spinor * const x,              /* (IN/OUT) initial guess on input, solution on output for this RHS*/
  spinor * const b               /* (IN) right-hand side*/
);

#endif 
