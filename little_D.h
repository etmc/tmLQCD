/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Recker, Carsten Urbach
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

#ifndef _LITTLE_D_H
#define _LITTLE_D_H

#include "complex.h"

extern int dfl_subspace_updated;
void little_D(complex * v, complex *w);
void unit_little_D(complex *v, complex *w);
void invert_little_D_spinor(spinor *r, spinor *s);
void apply_little_D_spinor(spinor *r, spinor *s);

#endif
