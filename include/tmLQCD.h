/***********************************************************************
 *
 * Copyright (C) 2014 Carsten Urbach
 *               2018 Bartosz Kostrzewa
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
 * invert wrapper for using tmLQCD as a library
 *
 * Author: Carsten Urbach
 *         curbach@gmx.de
 *
 *******************************************************************************/

#ifndef _TMLQCD_H
#define _TMLQCD_H

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct {
  unsigned int LX, LY, LZ, T, nstore, nsave, no_operators;
} tmLQCD_lat_params;

typedef struct {
  unsigned int nproc, nproc_t, nproc_x, nproc_y, nproc_z, cart_id, proc_id, time_rank,
      omp_num_threads;
  unsigned int proc_coords[4];
#ifdef TM_USE_MPI
    MPI_Comm cart_grid;
#else
    int cart_grid;
#endif
} tmLQCD_mpi_params;

typedef struct {
  // for the boundary phases
  // a value of 1.0 implies a theta angle of \pi/L
  double theta_x;
  double theta_y;
  double theta_z;
  double theta_t;

  // the twisted mass parameters are specified WITHOUT the 2kappa factor
  double mu;
  double mubar;
  double epsbar;
  double kappa;
  double c_sw;
} tmLQCD_op_params;

// in order to be able to initialise QMP if QPhiX is used, we need
// to allow tmLQCD to intialise MPI via QMP
// the input file has to be read because we set the number of threads as well
int tmLQCD_init_parallel_and_read_input(int argc, char* argv[], const int _verbose,
                                        char const * const input_filename);

int tmLQCD_invert_init(int argc, char *argv[], const int verbose, const int external_id);
int tmLQCD_read_gauge(const int nconfig);

// invert with source and propagator provided in TXYZ spin colour complex lexicographic order
// propagator has kappa normalisation
int tmLQCD_invert(double *const propagator, double *const source, const int op_id,
                  const int write_prop);

// invert on odd part of lattice with prepared source
int tmLQCD_invert_eo(double *const Odd_out, double *const Odd_in, const int op_id);

int tmLQCD_finalise();

int tmLQCD_get_gauge_field_pointer(double **gf);
int tmLQCD_get_mpi_params(tmLQCD_mpi_params *params);
int tmLQCD_get_lat_params(tmLQCD_lat_params *params);

int tmLQCD_set_op_params(tmLQCD_op_params const *const params, const int op_id);
int tmLQCD_get_op_params(tmLQCD_op_params *params, const int op_id);

  // need interfaces with (double*) pointers for these since we don't want to export
  // all of tmLQCDs data types
//// generates a point source at the global coordinates passed via the
//// four element vector global_txyz_src_pos in the ordering {t,x,y,z}
//// the spin index 'is' and the colour index 'ic'
//void full_source_spinor_field_point(spinor * const full_spinor,
//                                    const int is, const int ic,
//                                    const int * const global_txyz_src_pos);
//
//// as full_source_spinor_field_point but with output directly to checkerboarded
//// spinors
//void eo_source_spinor_field_point(spinor * const even_cb_spinor,
//                                  spinor * const odd_cb_spinor,
//                                  const int is, const int ic,
//                                  const int * const global_txyz_src_pos);
//
//void full_source_spinor_field_spin_diluted_oet_ts(spinor * const full_spinor,
//                                                  const int src_ts,
//                                                  const int src_d,
//                                                  const int sample,
//                                                  const int nstore,
//                                                  const unsigned int oet_seed);
//
//void eo_source_spinor_fied_spin_diluted_oet_ts(spinor * const even_cb_spinor,
//                                               spinor * const odd_cb_spinor,
//                                               const int src_ts,
//                                               const int src_d,
//                                               const int sample,   
//                                               const int nstore,                                            
//                                               const unsigned int oet_seed);

#ifdef TM_USE_QUDA
// FIXME: this needs to reworked into tmLQCD_invert, such that there is just
// one inversion pathway
// direct line to QUDA inverter, no messing about with even/odd reordering
// source and propagator  Should be full VOLUME spinor fields 
// op_id                  Index of the operator to be inverted (0 to N-1)
int invert_quda_direct(double *const propgator, double *const source, const int op_id);
#endif

#ifdef TM_USE_QPHIX
int tmLQCD_invert_qphix_direct(double *const Odd_out, double *const Odd_in, const int op_id);
#endif

#ifdef __cplusplus
}
#endif

#endif
