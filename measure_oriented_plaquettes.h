/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch, 2012 Bartosz Kostrzewa
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

#ifndef _MEASURE_ORIENTED_PLAQUETTES_H
#define _MEASURE_ORIENTED_PLAQUETTES_H

#include "su3.h"

/* measures the lattice average of plaquettes oriented in the 6
   hyperplanes TX, TY, TZ, XY, XZ, YZ and stores them in this
   order in the plaq array (of 6 elements) 
   
   the caller must provide the memory for plaq */

void measure_oriented_plaquettes(const su3 ** const gf, double *plaq);

/* implements the online measurement function for the oriented
   plaquettes, writes (in append mode) into "oriented_plaquettes.data" */

void oriented_plaquettes_measurement(const int traj, const int id, const int ieo);

#endif  
