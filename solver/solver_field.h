/***********************************************************************
 *
 * Copyright (C) 2009 Carsten Urbach
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

#ifndef _SOLVER_FIELD_H
#define _SOLVER_FIELD_H

#include "su3.h"

#include <complex.h>

int init_solver_field(spinor *** const solver_field, const int V, const int nr);
void finalize_solver(spinor ** solver_field, const int nr);
int init_solver_field_32(spinor32 *** const solver_field, const int V, const int nr);
void finalize_solver_32(spinor32 ** solver_field, const int nr);
int init_bisolver_field(bispinor *** const solver_field, const int V, const int nr);
void finalize_bisolver(bispinor ** solver_field, const int nr);
/* little solver fields (deflation) */
int init_lsolver_field(_Complex double *** const solver_field, const int V, const int nr);
void finalize_lsolver(_Complex double ** solver_field, const int nr);
int init_lsolver_field_32(_Complex float *** const solver_field, const int V, const int nr);
void finalize_lsolver_32(_Complex float ** solver_field, const int nr);

#endif
