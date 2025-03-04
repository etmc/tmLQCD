/***********************************************************************
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
#ifndef _EIGENVALUESJ_H
#define _EIGENVALUESJ_H

#include "su3.h"

extern su3_vector *eigenvectors_su3v;
extern double *eigenvls_su3v;
extern double * inv_eigenvls_su3v;
extern int eigenvalues_for_cg_computed_su3v;
extern int no_eigenvalues_su3v;
extern int evlength_su3v;

double eigenvalues_Jacobi(int * nr_of_eigenvalues, const int max_iterations, 
			  const double precision, const int maxmin, int tslice, const int nstore);

#endif // _EIGENVALUESJ_H
