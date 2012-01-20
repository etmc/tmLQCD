/***********************************************************************
 * Copyright (C) 2006 Thomas Chiarappa
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

#ifndef _ASSIGN_ADD_MUL_ADD_MUL_BI_H
#define _ASSIGN_ADD_MUL_ADD_MUL_BI_H

#include "su3.h"

/* (*R) = (*R) + c1*(*S) + c2*(*U) */
void assign_add_mul_add_mul_bi(bispinor * const R, bispinor * const S, bispinor * const U, const complex c1, const complex c2, const int N);


#endif
