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
#ifndef _POLYAKOV_LOOP_H
#define _POLYAKOV_LOOP_H

#include "measurements.h"

void polyakov_loop(_Complex double * pl_, const int mu);
int polyakov_loop_0(const int nstore, _Complex double* pl);
int polyakov_loop_dir(const int nstore, const int dir);
void polyakov_loop_measurement(const int nstore, const int id, const int ieo);

#endif
