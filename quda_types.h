/***********************************************************************
 * Copyright (C) 2017 Bartosz Kostrzewa
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
 *
 *******************************************************************************/

#ifndef TM_QUDA_TYPES_H
#define TM_QUDA_TYPES_H

#ifdef TM_USE_QUDA
#include <quda.h>
#else
// some definitions which are found in quda.h must be reproduced in case
// we are compiling without it, such that tm_QudaParams_t can be
// defined properly anyway
#define QUDA_MAX_MG_LEVEL 4
#include "quda_solver_translate.h"
#endif

typedef enum tm_quda_ferm_bc_t {
  TM_QUDA_THETABC = 0,
  TM_QUDA_APBC,
  TM_QUDA_PBC
} tm_quda_ferm_bc_t;


/* tm_QudaParams_t provides an interface between the tmLQCD input file and the
 * available QUDA parameters. At the moment, only the fermionic bounday conditions
 * and the MG parameters are exposed like this, but a further refactoring might
 * turn this into a complete representation of the possible input parameters */
typedef struct tm_QudaParams_t {
  tm_quda_ferm_bc_t fermionbc;

  int               mg_n_level;
  int               mg_n_vec[QUDA_MAX_MG_LEVEL];
  // in principle we can have mu scaling factors for all levels, but in practice
  // only the coarsest one is scaled
  double            mg_mu_factor;
  QudaInverterType    mg_setup_inv_type;
  int               mg_setup_maxiter;
} tm_QudaParams_t;

#endif // TM_QUDA_TYPES_H
