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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.	If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************/
/***********************************************************************
*
* File quda_interface.h
*
* Author: Mario Schroeck <mario.schroeck@roma3.infn.it>
* 
* Last changes: 04/2015
*
*
* Interface to QUDA for multi-GPU inverters
*
* The externally accessible functions are
*
*	 void _initQuda()
*		 Initializes the QUDA library. Carries over the lattice size and the
*		 MPI process grid and thus must be called after initializing MPI.
*		 Currently it is called in init_operators() if optr->use_qudainverter
*		 flag is set.
*		 Memory for the QUDA gaugefield on the host is allocated but not filled
*		 yet (the latter is done in _loadGaugeQuda(), see below).
*		 Performance critical settings are done here and can be changed.
*
*	 void _endQuda()
*		 Finalizes the QUDA library. Call before MPI_Finalize().
*
*	 void _loadGaugeQuda()
*		 Copies and reorders the gaugefield on the host and copies it to the GPU.
*		 Must be called between last changes on the gaugefield (smearing etc.)
*		 and first call of the inverter. In particular, 'boundary(const double kappa)'
*		 must be called before if nontrivial boundary conditions are to be used since
*		 those will be applied directly to the gaugefield. Currently it is called just
*		 before the inversion is done (might result in wasted loads...).
*
*	 double tmcgne_quda(int nmx,double res,int k,int l,int *status,int *ifail)
*		 The same functionality as 'tmcgne' (see tmcg.c) but inversion is performed on
*		 the GPU using QUDA. Final residuum check is performed on the host (CPU)
*		 with the function 'void tmQnohat_dble(int k,int l)' (see tmdirac.c).
*
*	 The functions
*
*		 int invert_eo_quda(...);
*		 int invert_clover_eo_quda(...);
*		 int invert_doublet_eo_quda(...);
*		 int invert_cloverdoublet_eo_quda(...);
*		 void M_full_quda(...);
*		 void D_psi_quda(...);
*
*	 closely mimic their tmLQCD counterparts in functionality as well as
*	 input and output parameters.
*
*	 To activate those, set "UseQudaInverter = yes" in the operator
*	 declaration of the input file. For details see the documentation.
*
*
* Notes:
*
* Minimum QUDA version is 0.7.0 (see https://github.com/lattice/quda/issues/151 
* and https://github.com/lattice/quda/issues/157).
*
*
**************************************************************************/

#ifndef QUDA_INTERFACE_H_
#define QUDA_INTERFACE_H_

#include "su3.h"
#include "solver/solver_params.h"

// pointer to the QUDA gaugefield
double *gauge_quda[4];

// wrapper functions
void _initQuda();
void _endQuda();
void _loadGaugeQuda();

// to be called instead of tmcgne to use the QUDA inverter
int invert_eo_quda(spinor * const Even_new, spinor * const Odd_new,
                   spinor * const Even, spinor * const Odd,
                   const double precision, const int max_iter,
                   const int solver_flag, const int rel_prec,
                   const int sub_evs_flag, const int even_odd_flag,
                   const int no_extra_masses, double * const extra_masses, solver_params_t solver_params,
                   const int id );

int invert_clover_eo_quda(spinor * const Even_new, spinor * const Odd_new,
                          spinor * const Even, spinor * const Odd,
                          const double precision, const int max_iter,
                          const int solver_flag, const int rel_prec, solver_params_t solver_params);

int invert_doublet_eo_quda(spinor * const Even_new_s, spinor * const Odd_new_s,
                           spinor * const Even_new_c, spinor * const Odd_new_c,
                           spinor * const Even_s, spinor * const Odd_s,
                           spinor * const Even_c, spinor * const Odd_c,
                           const double precision, const int max_iter,
                           const int solver_flag, const int rel_prec);

int invert_cloverdoublet_eo_quda(spinor * const Even_new_s, spinor * const Odd_new_s,
                           spinor * const Even_new_c, spinor * const Odd_new_c,
                           spinor * const Even_s, spinor * const Odd_s,
                           spinor * const Even_c, spinor * const Odd_c,
                           const double precision, const int max_iter,
                           const int solver_flag, const int rel_prec);

int invert_quda(spinor * const P, spinor * const Q, const int max_iter, double eps_sq, const int rel_prec );

// apply the TM operator using QUDA
void M_full_quda(spinor * const Even_new, spinor * const Odd_new,  spinor * const Even, spinor * const Odd);
void D_psi_quda(spinor * const P, spinor * const Q);

#endif /* QUDA_INTERFACE_H_ */
