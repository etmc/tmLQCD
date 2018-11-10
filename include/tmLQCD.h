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

#ifdef TM_USE_QUDA
  int invert_quda_direct(double * const propgator, double * const source,
                    const int op_id);
#endif

#ifdef __cplusplus
}
#endif

#endif
