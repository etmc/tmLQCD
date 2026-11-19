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
#ifndef _MEASURE_RECTANGLES_H
#define _MEASURE_RECTANGLES_H

#include "su3.h"

/* apply_ptbc == 1 weights each rectangle with the product of the PTBC defect
 * coefficients of its six links (use for the gauge action), 0 gives the plain
 * rectangle observable. */
double measure_rectangles(const su3** const gf, int const apply_ptbc);

#endif
