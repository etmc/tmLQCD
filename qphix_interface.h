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

#include "su3.h"

typedef enum QphixPrec { QPHIX_FLOAT_PREC = 0, QPHIX_HALF_PREC, QPHIX_DOUBLE_PREC } QphixPrec;

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

// Initialize and Finalize QPhiX
void _initQphix(int argc, char **argv, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_,
                int PadXYZ_, int MinCt_, int c12, QphixPrec precision_);
void _endQphix();

// Wrapper functions for Full Solver and Dslash
void invert_qphix(spinor *const P, spinor *const Q, const int max_iter, double eps_sq,
                 const int rel_prec);
void D_psi_qphix(spinor* Odd_out, const spinor* Odd_in);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif
#endif /* QPHIX_INTERFACE_H_ */
