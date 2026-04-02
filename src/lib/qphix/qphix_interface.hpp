/***********************************************************************
 *
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
 ***********************************************************************/

#pragma once

#include "global.h"
#include "qphix_types.h"

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

#include "misc_types.h"
#include "operator_types.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/solver_params.h"
#include "su3.h"

#ifdef __cplusplus
}
#endif

#include <vector>

int invert_eo_qphix_nflavour_mshift(std::vector< std::vector< spinor* > > &Odd_out, 
                                    std::vector< std::vector< spinor* > > &Odd_in, 
                                    const double precision,
                                    const int max_iter,
                                    const int solver_flag, 
                                    const int rel_prec,
                                    solver_params_t solver_params,
                                    const SloppyPrecision sloppy, const CompressionType compression,
                                    const int num_flavour);