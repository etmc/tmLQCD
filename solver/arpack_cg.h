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

#ifndef _ARPACK_CG_H
#define _ARPACK_CG_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_params.h"
#include "solver/eigenvalues_arpack.h"


int arpack_cg(
  /* solver params */
  const int N,                   /* (IN) Number of lattice sites for this process*/
  solver_params_t solver_params, /* (IN) parameters for solver */
  spinor * const x,              /* (IN/OUT) initial guess on input, solution on output for this RHS*/
  spinor * const b,              /* (IN) right-hand side*/
  matrix_mult f,                 /* (IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply in double precision */
  matrix_mult f32,               /* (IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply in single precision */
  const double eps_sq,           /* (IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
  const int rel_prec,            /* (IN) 0 for using absoute error for convergence
                                         1 for using relative error for convergence*/
  const int maxit,               /* (IN) Maximum allowed number of iterations to solution for the linear system*/
  matrix_mult f_final,           /* (IN) final operator application during projection of type 1 */
  matrix_mult f_initial          /* (IN) initial operator application during projection of type 1 */
);

#endif 
