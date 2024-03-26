/***********************************************************************
 *
 * Copyright (C) 2024 Bartosz Kostrzewa
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
 ************************************************************************/

#ifdef HAVE_CONFIG_H
# include <tmlqcd_config.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
#endif

#include <string.h>
#include <stdio.h>

#include "global.h"
#include "geometry_eo.h"
#include "fatal_error.h"
#include "measurements.h"
#include "operator.h"
#include "gettime.h"
#include "tm_debug_printf.h"
#include "misc_types.h"

#ifdef TM_USE_QUDA
#include "quda_interface.h"
#endif

void eigenvalues_measurement(const int traj, const int id, const int ieo) {
  tm_stopwatch_push(&g_timers, __func__, "");
  init_operators();

  if(no_operators < 1){
    tm_debug_printf(0, 0, "Error: no operators defined in input file, cannot perform eigenvalues online measurement!\n");
#ifdef TM_USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  eig_param_t eig = measurement_list[id].eig;

  double * evals = calloc(eig.n_evals, sizeof(double));

  for( int op_id = 0; op_id < no_operators; op_id++ ){
    operator * optr = &operator_list[op_id];

    if( !(optr->type == TMWILSON || optr->type == WILSON ||
          optr->type == CLOVER || optr->type == DBTMWILSON ||
          optr->type == DBCLOVER) ){
      tm_debug_printf(0, 0, "Error: only operator types WILSON, TMWILSON, CLOVER, DBTMWILSON and DBCLOVER are supported.\n");
#ifdef TM_USE_MPI
      MPI_Finalize();
#endif
      exit(1);
    }

    op_backup_restore_globals(TM_BACKUP_GLOBALS);
    op_set_globals(op_id);

    if( measurement_list[id].external_library == QUDA_LIB ){
#ifdef TM_USE_QUDA
      eigsolveQuda(evals, eig.n_evals, eig.tol, 1, 0, eig.max_iter, eig.maxmin,
                   optr->eps_sq, optr->maxiter, eig.polydeg, eig.amin, eig.amax, eig.n_kr,
                   optr->solver, optr->solver, optr->rel_prec, ieo,
                   optr->sloppy_precision, optr->sloppy_precision, optr->compression_type,
                   optr->type == CLOVER || optr->type == TMWILSON || optr->type == WILSON); 
#else
    tm_debug_printf(0, 0, "Error: Attempted to use QUDA eigensolver but this build was not configured for QUDA usage.\n");
#ifdef TM_USE_MPI
    MPI_Finalize();
#endif
    exit(1);
#endif
    }
    
    op_backup_restore_globals(TM_RESTORE_GLOBALS);
  
  } // loop over operators
  
  free(evals);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}

