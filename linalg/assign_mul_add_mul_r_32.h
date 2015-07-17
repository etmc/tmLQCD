/***********************************************************************
 * Copyright (C) 2015 Florian Burger
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

#ifndef _ASSIGN_MUL_ADD_MUL_R_32_H
#define _ASSIGN_MUL_ADD_MUL_R_32_H

#include "su3.h"

/* Makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
void assign_mul_add_mul_r_32(spinor32 * const R,spinor32 * const S, 
			  const float c1, const float c2,
			  const int N);

#endif
