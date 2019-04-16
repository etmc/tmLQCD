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

#ifndef QUDA_DUMMY_TYPES_H
#define QUDA_DUMMY_TYPES_H

#include "solver/solver_types.h"

// some definitions which are found in quda.h must be reproduced in case
// we are compiling without it, such that tm_QudaParams_t can be
// defined properly anyway
#define QUDA_MAX_MG_LEVEL 4
#define QUDA_BOOLEAN_YES 1
#define QUDA_BOOLEAN_NO 0

#include <limits.h>
#define QUDA_INVALID_ENUM INT_MIN

// these exist only in case we are compiling without QUDA support, such that the
// input file reader can be compiled
typedef enum QudaInverterType_s {
 QUDA_BICGSTAB_INVERTER = BICGSTAB,
 QUDA_CG_INVERTER = CG,
 QUDA_MR_INVERTER = MR,
 QUDA_GCR_INVERTER = GCR,
 QUDA_CA_GCR_INVERTER = CA_GCR
} QudaInverterType;

typedef QudaCABasis_s {
  QUDA_POWER_BASIS,
  QUDA_CHEBYSHEV_BASIS,
  QUDA_INVALID_BASIS = QUDA_INVALID_ENUM
} QudaCABasis;

#endif
