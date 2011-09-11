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
#ifndef _CG_HERSU3V_H
#define _CG_HERSU3V_H

#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int cg_her_su3vect(su3_vector * const P, su3_vector * const Q, const int max_iter, double eps_sq, const int rel_prec, 
		   const int N, const int tslice, matrix_mult_su3vect f);

#endif
