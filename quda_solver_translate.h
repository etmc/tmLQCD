/***************************************************************************
 * Copyright (C) 2017                               Bartosz Kostrzewa
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
 ****************************************************************************/

#ifndef QUDA_SOLVER_TRANSLATE_H
#define QUDA_SOLVER_TRANSLATE_H

#include "solver/solver_types.h"

// these exist only in case we are compiling without QUDA support, such that the
// input file reader can be compiled
typedef enum QudaInverterType_s {
 QUDA_BICGSTAB_INVERTER = BICGSTAB,
 QUDA_CG_INVERTER = CG
} QudaInverterType;

#endif
