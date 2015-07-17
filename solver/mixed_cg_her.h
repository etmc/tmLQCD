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
#ifndef _MIXED_CG_HER_H
#define _MIXED_CG_HER_H

#include"operator/tm_operators_32.h"
#include"solver/matrix_mult_typedef.h"
#include"su3.h"

int mixed_cg_her(spinor * const P, spinor * const Q, const int max_iter, 
		 double eps_sq, const int rel_prec, const int N, matrix_mult f, matrix_mult32 f32);

#endif
