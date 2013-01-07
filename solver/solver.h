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

#ifndef _SOLVER_H
#define _SOLVER_H

#define BICGSTAB 0
#define CG 1
#define GMRES 2
#define CGS 3
#define MR 4
#define BICGSTABELL 5
#define FGMRES 6
#define GCR 7
#define GMRESDR 8
#define PCG 9
#define DFLGCR 10
#define DFLFGMRES 11
#define CGMMS 12
#define MIXEDCG 13
#define CGMMSND 14

#include"solver/matrix_mult_typedef.h"

#include"solver/gmres.h"
#include"solver/gmres_dr.h"
#include"solver/fgmres.h"
#include"solver/bicgstab_complex.h"
#include"solver/cgs_real.h"
#include"solver/bicgstabell.h"
#include"solver/bicgstab2.h"
#include"solver/cg_her.h"
#include"solver/pcg_her.h"
#include"solver/mr.h"
#include"solver/gcr.h"
#include"solver/eigenvalues.h"
#include"solver/cg_mms_tm.h"
#include"solver/mixed_cg_her.h"

#include"solver/sub_low_ev.h"
#include"solver/gmres_precon.h"
#include"solver/poly_precon.h"

#include "solver/matrix_mult_typedef_bi.h"
#include "solver/bicgstab_complex_bi.h"
#include "solver/cg_her_bi.h"

#include "solver/matrix_mult_typedef_nd.h"
#include "solver/cg_her_nd.h"

#include "solver/generate_dfl_subspace.h"

typedef struct {
  int max_iter, rel_prec, type, no_shifts, N;
  double eps_sq;
  matrix_mult f;
  matrix_mult_nd g;
  double * shifts;
} solver_pm_t;
#endif
