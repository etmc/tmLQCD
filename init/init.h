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

#ifndef _INIT_H
#define _INIT_H

#include "init/init_bispinor_field.h"
#include "init/init_chi_spinor_field.h"
#include "init/init_dirac_halfspinor.h"
#include "init/init_gauge_field.h"
#include "init/init_gauge_tmp.h"
#include "init/init_geometry_indices.h"
#ifdef WITHLAP
#  include "init/init_jacobi_field.h"
#endif
#include "init/init_moment_field.h"
#include "init/init_spinor_field.h"
#include "init/init_stout_smear_vars.h"
#ifdef OMP
# include <omp.h>
# include "init/init_omp_accumulators.h"
# include "init/init_openmp.h"
#endif

#endif
