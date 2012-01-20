/***********************************************************************
 *
 * Copyright (C) 2005 Martin Hasenbusch
 *               2009 Carsten Urbach
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

#ifndef _CLOVER_H
#define _CLOVER_H

#include "su3.h"

extern su3 *** sw;
extern su3 *** sw_inv;
extern su3 ** swm, ** swp;
extern su3adj ** dclover;

void clover_inv(const int ieo, spinor * const l);
void clover_gamma5(const int ieo, 
		   spinor * const l, spinor * const k, spinor * const j,
		   const double q_off);
void clover(const int ieo, 
	    spinor * const l, spinor * const k, spinor * const j,
	    const double q_off);

#endif
