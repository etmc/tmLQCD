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

#ifndef _CHRONO_GUESS_H
#define _CHRONO_GUESS_H

#include "solver/matrix_mult_typedef.h"

void chrono_add_solution(spinor * const trial, spinor ** const v, int index_array[],
			const int _N, int * _n, const int V);

int chrono_guess(spinor * const trial, spinor * const phi, spinor ** const v, int index_array[], 
		 const int N, const int n, const int V, matrix_mult f);

#endif
