/***************************************************************************
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
 ****************************************************************************/

/*****************************************************************************
 * Struct for passing the parameters for the solver. This should replace all 
 * solver related parameters and eliminate the need for global parameters as 
 * this struct will be a member of the operator struct.
 * 
 * A. M. Abdel-Rehim (a.abdel-rehim@cyi.ac.cy)
 * March, 17th, 2013
 ****************************************************************************/


#ifndef _SOLVER_PARAMS_H
#define _SOLVER_PARAMS_H

typedef struct {

  /********************************
   * Incremental EigCG parameters
   ********************************/

  int eigcg_nrhs;          /*total number of right-hand sides to be solved*/
  int eigcg_nev;           /*number of eigenvalues computed from a single right-hand side */
  int eigcg_vmax;          /*size of the search subspace for eigcg*/
  int eigcg_ldh;           /*total number of eigenvectors that will be computed and used in deflation */
  double eigcg_restolsq;   /*tolerance squared for restarting eigcg after eigenvectors has been computed
                             Typically this is the square root of the tolerance squared requested for the linear system.
                             Example, to solve the linear systems to res_sq=1e-16, one chooses eigcg_restolsq=1e-8 or smaller 
                             This will specify how many times deflated CG restaretd in the second phase (after eigenvectors has been computed)*/
} solver_params_t;



#endif


 
