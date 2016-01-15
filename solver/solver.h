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


#include"solver/solver_types.h"
#include"solver/matrix_mult_typedef.h"
#include "solver/matrix_mult_typedef_bi.h"
#include "solver/matrix_mult_typedef_nd.h"

typedef struct {
  // solver type
  int type;
  // maximal number of iterations
  int max_iter;
  // use relative precision
  int rel_prec;
  // number of shifts in multi shift solvers
  int no_shifts;
  // dimension of spinors
  int sdim;
  // squared desired residue
  double squared_solver_prec;
  // single flavour matrix to invert
  matrix_mult M_psi;
  // 32bit single flavour matrix to invert
  matrix_mult32 M_psi32;  
  // flavour doublet matrix to invert
  matrix_mult_nd M_ndpsi;
  // 32bit flavour doublet matrix to invert
  matrix_mult_nd32 M_ndpsi32;  
  // pointer to array of shifts
  double * shifts;
} solver_pm_t;

#include"solver/gmres.h"
#include"solver/gmres_dr.h"
#include"solver/fgmres.h"
#include"solver/bicgstab_complex.h"
#include"solver/bicg_complex.h"
#include"solver/cgs_real.h"
#include"solver/bicgstabell.h"
#include"solver/bicgstab2.h"
#include"solver/cg_her.h"
#include"solver/pcg_her.h"
#include"solver/mr.h"
#include"solver/gcr.h"
#include"solver/incr_eigcg.h"
#include"solver/eigenvalues.h"
#include"solver/cg_mms_tm.h"
#include"solver/mixed_cg_her.h"
#include"solver/mcr.h"
#include"solver/cr.h"

#include"solver/sub_low_ev.h"
#include"solver/gmres_precon.h"
#include"solver/poly_precon.h"

#include "solver/bicgstab_complex_bi.h"
#include "solver/cg_her_bi.h"

#include "solver/cg_her_nd.h"
#include"solver/cg_mms_tm_nd.h"
#include"solver/mixed_cg_mms_tm_nd.h"

#include "solver/generate_dfl_subspace.h"

#include "solver/sumr.h"

#endif
