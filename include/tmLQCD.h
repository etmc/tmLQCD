/***********************************************************************
 *
 * Copyright (C) 2014 Carsten Urbach
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
#include "su3.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
  
  typedef struct {
    unsigned int LX, LY, LZ, T, nstore, nsave, no_operators;
  } tmLQCD_lat_params;

  typedef struct {
    unsigned int nproc, nproc_t, nproc_x, nproc_y, nproc_z, cart_id, proc_id, time_rank, omp_num_threads;
    unsigned int proc_coords[4];
  } tmLQCD_mpi_params;

  int tmLQCD_invert_init(int argc, char *argv[], const int verbose, const int external_id);
  int tmLQCD_read_gauge(const int nconfig);
  int tmLQCD_invert(double * const propagator, double * const source,
		    const int op_id, const int write_prop);
  int tmLQCD_finalise();

  int tmLQCD_get_gauge_field_pointer(double ** gf);
  int tmLQCD_get_mpi_params(tmLQCD_mpi_params * params);
  int tmLQCD_get_lat_params(tmLQCD_lat_params * params);

// generates a point source at the global coordinates passed via the
// four element vector global_txyz_src_pos in the ordering {t,x,y,z}
// the spin index 'is' and the colour index 'ic'
void full_source_spinor_field_point(spinor * const full_spinor,
                                    const int is, const int ic,
                                    const int * const global_txyz_src_pos);

// as full_source_spinor_field_point but with output directly to checkerboarded
// spinors
void eo_source_spinor_field_point(spinor * const even_cb_spinor,
                                  spinor * const odd_cb_spinor,
                                  const int is, const int ic,
                                  const int * const global_txyz_src_pos);

void full_source_spinor_field_spin_diluted_oet_ts(spinor * const full_spinor,
                                                  const int src_ts,
                                                  const int src_d,
                                                  const int sample,
                                                  const int nstore,
                                                  const unsigned int oet_seed);

void eo_source_spinor_fied_spin_diluted_oet_ts(spinor * const even_cb_spinor,
                                               spinor * const odd_cb_spinor,
                                               const int src_ts,
                                               const int src_d,
                                               const int sample,   
                                               const int nstore,                                            
                                               const unsigned int oet_seed);

#ifdef TM_USE_QUDA
  int invert_quda_direct(double * const propgator, double * const source,
                    const int op_id);
#endif

#ifdef __cplusplus
}
#endif

#endif
