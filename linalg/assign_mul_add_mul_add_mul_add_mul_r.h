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

#ifndef _ASSIGN_MUL_ADD_MUL_ADD_MUL_ADD_MUL_R_H
#define _ASSIGN_MUL_ADD_MUL_ADD_MUL_ADD_MUL_R_H

#include "su3.h"

/* Makes (*R) = c1*(*R) + c2*(*S) + c3*(*U) + c4*(*V)*/
void assign_mul_add_mul_add_mul_add_mul_r(spinor * const R,spinor * const S,spinor * const U,spinor * const V,const double c1,const double c2,const double c3,const double c4, const int N);


#endif
