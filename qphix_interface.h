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
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>
*
* Last changes: 03/2015
*
*
* Integration of the QUDA inverter for multi-GPU usage
*
* The externally accessible functions are
*
*   void _initQphix( int verbose )
*     Initializes the QUDA library. Carries over the lattice size and the
*     MPI process grid and thus must be called after initializing MPI (and
*     after 'read_infile(argc,argv)').
*     Memory for the QUDA gaugefield on the host is allocated but not filled
*     yet (the latter is done in _loadGaugeQphix(), see below).
*     Performance critical settings are done here and can be changed.
*     Input parameter: verbose (0=SILENT, 1=SUMMARIZE, 2=VERBOSE).
*
*   void _endQphix()
*     Finalizes the QUDA library. Call before MPI_Finalize().
*
*   void _loadGaugeQphix()
*     Copies and reorders the gaugefield on the host and copies it to the GPU.
*     Must be called between last changes on the gaugefield (smearing, flip
*     boundary conditions, etc.) and first call of the inverter.
*
*   double tmcgne_quda(int nmx,double res,int k,int l,int *status,int *ifail)
*     The same functionality as 'tmcgne' (see tmcg.c) but inversion is performed
*on
*     the GPU using QUDA. Final residuum check is performed on the host (CPU)
*     with the function 'void tmQnohat_dble(int k,int l)' (see tmdirac.c).
*
*   void tmQnohat_quda(int k, int l)
*     The implementation of the QUDA equivalent of 'tmQnohat_dble'.
*
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/quda/issues/151
* and https://github.com/lattice/quda/issues/157).
*
* To enable compilation of the same code for QUDA usage and standard non-QUDA
*usage,
* all calls of these functions should be wrapped in precompiler switches of the
*form
*
*   #ifdef QUDA
*     ...
*   #endif
*
**************************************************************************/

#ifndef QPHIX_INTERFACE_H_
#define QPHIX_INTERFACE_H_

#include "su3.h"

double *gauge_qphix[4];

typedef enum QphixPrec { QPHIX_FLOAT_PREC = 0, QPHIX_HALF_PREC, QPHIX_DOUBLE_PREC } QphixPrec;

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

// wrapper functions
void _initQphix(int argc, char **argv, int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_,
                int PadXYZ_, int MinCt_, int c12, QphixPrec precision_);
void _endQphix();
void _loadGaugeQphix();

// to be called instead of tmcgne to use the QUDA inverter
int invert_qphix(spinor *const P, spinor *const Q, const int max_iter, double eps_sq,
                 const int rel_prec);

// apply the TM operator using QUDA
void M_full_qphix(spinor *const Even_new, spinor *const Odd_new, spinor *const Even,
                  spinor *const Odd);
void D_psi_qphix(spinor* Odd_out, const spinor* Odd_in);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif
#endif /* QPHIX_INTERFACE_H_ */
