/*****************************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim. The original code was written
 * by Andreas Stathopoulos and Kostas Orginos and integrated in Chroma.
 * In this interface we use functions from tmLQCD.
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
 *
 * Incremental eigCG for solving linear systems multiple right-hand sides
 ****************************************************************************/
/* A sample input is given in the sample-input folder */

#ifndef _INCR_EIGCG_H
#define _INCR_EIGCG_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"


int incr_eigcg(
     const int N,      /*(IN) Number of lattice sites for this process*/
     const int nrhs,   /*(IN) Number of right-hand sides to be solved*/ 
     const int nrhs1,  /*(IN) First number of right-hand sides to be solved using tolerance eps_sq1*/ 
     spinor * const x, /*(IN/OUT) initial guess on input, solution on output for this RHS*/
     spinor * const b, /*(IN) right-hand side*/
     const int ldh,    /*(IN) maximum number of eignvectors to be computed*/
     matrix_mult f,    /*(IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply*/
     const double eps_sq1,    /*(IN) squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
     const double eps_sq,     /*(IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
     double restart_eps_sq, /*(IN) squared tolerance for restarting CG*/
     const int rand_guess_opt,    /*(IN) set to non-zero if you want to use random intitial guess (volume Gaussian with mean 0)*/
     const int rel_preq,    /*(IN)0 for using absoute error for convergence
                                  1 for using relative error for convergence*/
     const int maxit,       /*(IN) Maximum allowed number of iterations to solution*/
     int nev,               /*(IN)number of eigenvectors to be computed while solving
                                  every right-hand side until the maximum number ldh is reached*/
     const int v_max);  /*(IN) subspace size used to compute nev vectors*/

#endif
