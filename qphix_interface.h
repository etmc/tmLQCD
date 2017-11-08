/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
 *               2017 Peter Labus, Martin Ueding, Bartosz Kostrzewa
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

#ifndef QPHIX_INTERFACE_H_
#define QPHIX_INTERFACE_H_

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

int invert_eo_qphix_oneflavour(spinor* const Odd_out, spinor* const Odd_in, const int max_iter,
                               const double precision, const int solver_flag, const int rel_prec,
                               solver_params_t solver_params, const SloppyPrecision sloppy,
                               const CompressionType compression);

int invert_eo_qphix_oneflavour_mshift(spinor** Odd_out, spinor* const Odd_in, const int max_iter,
                                      const double precision, const int solver_flag, const int rel_prec,
                                      solver_params_t solver_params, const SloppyPrecision sloppy,
                                      const CompressionType compression);

int invert_eo_qphix_twoflavour(spinor* Odd_out_s, spinor* Odd_out_c, spinor* Odd_in_s,
                               spinor* Odd_in_c, const int max_iter, const double precision,
                               const int solver_flag, const int rel_prec,
                               solver_params_t solver_params, const SloppyPrecision sloppy,
                               const CompressionType compression);

int invert_eo_qphix_twoflavour_mshift(spinor** Odd_out_s, spinor** Odd_out_c, spinor* Odd_in_s,
                                      spinor* Odd_in_c, const int max_iter, const double precision,
                                      const int solver_flag, const int rel_prec,
                                      solver_params_t solver_params, const SloppyPrecision sloppy,
                                      const CompressionType compression);

void Mfull_qphix(spinor* Even_out, spinor* Odd_out, const spinor* Even_in, const spinor* Odd_in,
                 const op_type_t op_type);

void testSpinorPackers(spinor* Even_out, spinor* Odd_out, const spinor* const Even_in,
                       const spinor* const Odd_in);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif
#endif /* QPHIX_INTERFACE_H_ */
