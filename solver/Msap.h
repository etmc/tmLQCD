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

#ifndef _MSAP_H
#define _MSAP_H

#include "su3.h"

void Msap(spinor * const P, spinor * const Q, const int Ncy, const int Niter);
void Msap_eo(spinor * const P, spinor * const Q, const int Ncy, const int Niter);
void CGeoSmoother(spinor * const P, spinor * const Q, const int Ncy, const int dummy);
void Mtm_plus_block_psi(spinor * const l, spinor * const k, const int i);
void Mtm_plus_sym_block_psi(spinor * const l, spinor * const k, const int i);
void Msw_plus_block_psi(spinor * l, spinor *  k, const int i);
//This function doesn't work properly yet because we need the inverse of 1+Too+img5 where
//Too is the clover term in the symmetric case
void Msw_plus_sym_block_psi(spinor *  l, spinor *  k, const int i);

#endif
