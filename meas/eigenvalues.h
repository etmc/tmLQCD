/***********************************************************************
 *
 * Copyright (C) 2024 Bartosz Kostrzewa
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

#ifndef _MEASURE_EIGENVALUES_H
#define _MEASURE_EIGENVALUES_H

#include "su3.h"

/* implements the eigenvalue online measurement
 * writes into eigenvalues.mXY.IJKLMN 
 * where XY corresponds to the operator index and IJKLMN refers to the trajectory counter */

void eigenvalues_measurement(const int traj, const int id, const int ieo);

#endif  
