/***********************************************************************
 * Copyright (C) 2015 Bartosz Kostrzewa
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
#ifndef _RG_MIXED_CG_HER_H
#define _RG_MIXED_CG_HER_H

#include"operator/tm_operators_32.h"
#include"solver/matrix_mult_typedef.h"
#include"solver/solver_params.h"
#include"su3.h"

typedef enum PolakRibiere_s {
  MCG_NO_PR=0,
  MCG_PR=1
} MCG_PR_TYPE;

typedef enum Pipelined_s {
  MCG_NO_PIPELINED=0,
  MCG_PIPELINED
} MCG_PIPELINED_TYPE;

// currently not used
typedef enum ResGuide_s {
  MCG_NO_RESGUIDE=0,
  MCG_RESGUIDE
} MCG_RESGUIDE_TYPE;

int rg_mixed_cg_her(spinor * const P, spinor * const Q, solver_params_t solver_params,
                    const int max_iter, const double eps_sq, const int rel_prec,
                    const int N, matrix_mult f, matrix_mult32 f32);

#endif
