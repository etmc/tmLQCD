/***********************************************************************
 * Copyright (C) 2010 Florian Burger
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


#ifndef _MIXED_SOLVE_H_

void initnn();
extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps,const int rel_prec, const int N);
extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N);

extern "C" int bind_texture_spin(dev_spinor* s, int i);
extern "C" int unbind_texture_spin(int i);


#define _MIXED_SOLVE_H_
#endif
