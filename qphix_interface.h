/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
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
/***********************************************************************
*
* File qphix_interface.h
*
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>,
*         Peter Labus <Peter.Labus@sissa.it>
*
* Last changes: 03/2017
*
*
* Integration of the QPhiX library for Intel Xeon Phi usage
*
* The externally accessible functions are
*
*   void _initQphix(int argc, char **argv,
*                   int By_, int Bz_, int NCores_,
*                   int Sy_, int Sz_, int PadXY_,
*                   int PadXYZ_, int MinCt_, int c12, QphixPrec precision_)
*     Initializes the QPhiX library. Carries over the lattice size and the
*     MPI process grid and thus must be called after initializing MPI (and
*     after 'read_infile(argc,argv)').
*
*   void _endQphix()
*     Finalizes the QPhiX library. Call before MPI_Finalize().
* Integration of the QUDA inverter for multi-GPU usage
*
*
* Notes:
*
* To enable compilation of the same code for QPhiX usage and standard non-QPhiX
* usage, all calls of these functions should be wrapped in precompiler switches
* of the form
*
*   #ifdef TM_USE_QPHIX
*     ...
*   #endif
*
**************************************************************************/

#ifndef QPHIX_INTERFACE_H_
#define QPHIX_INTERFACE_H_

#include "global.h"
#include "qphix_types.h"
#include "solver/solver_params.h"
#include "su3.h"
#include "operator_types.h"

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

// Initialize and Finalize QPhiX
void _initQphix(int argc, char** argv, QphixParams_t params, int c12, QphixPrec_t precision);
void _endQphix();

// Wrapper functions for Full Solver and Dslash
int invert_eo_qphix(spinor* const Even_new, spinor* const Odd_new, spinor* const Even,
                    spinor* const Odd, const double precision, const int max_iter,
                    const int solver_flag, const int rel_prec, solver_params_t solver_params,
                    const SloppyPrecision sloppy, const CompressionType compression);
void Mfull_qphix(spinor* Even_out, spinor* Odd_out, const spinor* Even_in, 
                 const spinor* Even_odd, const op_type_t op_type);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif
#endif /* QPHIX_INTERFACE_H_ */
