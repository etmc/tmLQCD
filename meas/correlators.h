/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
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

#ifndef _ONLINE_MEASUREMENT_H
#define _ONLINE_MEASUREMENT_H

#include <stdbool.h>

/* measurement of the correlators involving the light doublet (see tmLQCD documentation)*/
void light_correlators_measurement(const int traj, const int id, const int ieo);

/* measurement of the 4x4 matrix for the heavy doublet: eq. 20 of https://arxiv.org/pdf/1005.2042.pdf  */
void heavy_correlators_measurement(const int traj, const int id, const int ieo, const int i1, const int i2);

/* 
 Function that is called when the correlator measurement is specified in the input file.
 Internally, it calls the functions for the measure of the light and (optionally) the heavy mesons correlators
 */
void correlators_measurement(const int traj, const int id, const int ieo);


#endif
