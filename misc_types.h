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
 ***********************************************************************/

#ifndef MISC_TYPES_H
#define MISC_TYPES_H

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif

#ifdef TM_USE_MPI
#include <mpi.h>
#endif

#ifdef TM_USE_QPHIX
#include "qphix/qphix_config.h"
#endif

#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif

#include "tm_debug_printf.h"

#define TM_GAUGE_FIELD_NAME_LENGTH 100

// this number is completely arbitrary, but it's supposed to make sure that
// for example the QUDA MG setup is reset when a trajectory was rejected in the HMC 
#define TM_GAUGE_PROPAGATE_THRESHOLD 3.0
// this number is used to advance the state of the gauge field when SU3-restoration
// is done
#define TM_GAUGE_PROPAGATE_MIN 0.001

/* enumeration type for the identity of the program
 * which is being executed
 * this is useful to unify various utility functions which
 * otherwise lead to a lot of code duplication */
typedef enum tm_ProgramId_t {
  TM_PROGRAM_HMC_TM = 0,
  TM_PROGRAM_INVERT,
  TM_PROGRAM_OFFLINE_MEASUREMENT,
  TM_PROGRAM_BENCHMARK,
  TM_PROGRAM_DERIV_MG_TUNE,
  TM_PROGRAM_EXTERNAL
} tm_ProgramId_t;

/* enumeration type for return value 
 * we follow http://tldp.org/LDP/abs/html/exitcodes.html for the starting 
 * value */
typedef enum tm_ExitCode_t {
  TM_EXIT_SUCCESS = 0,
  TM_EXIT_INVALID_CMDLINE_ARG = 166
} tm_ExitCode_t;

/* enumeration type for the sloppy prec. of the inverter */
typedef enum SloppyPrecision_s {
  SLOPPY_DOUBLE = 0,
  SLOPPY_SINGLE,
  SLOPPY_HALF
} SloppyPrecision;

/* enumeration type for the compression of the inverter */
typedef enum CompressionType_s {
  NO_COMPRESSION = 18,
  COMPRESSION_12 = 12,
  COMPRESSION_8  = 8
} CompressionType;

/* enumeration type for the external inverter */
typedef enum ExternalInverter_s {
  NO_EXT_INV = 0,
  QUDA_INVERTER,
  QPHIX_INVERTER
} ExternalInverter;

/* enumeration type for the external inverter */
typedef enum ExternalLibrary_s {
  NO_EXT_LIB = 0,
  QUDA_LIB
} ExternalLibrary;


typedef enum backup_restore_t {
  TM_BACKUP_GLOBALS = 0,
  TM_RESTORE_GLOBALS
} backup_restore_t;

typedef enum real_imag_t {
  TM_REAL = 0,
  TM_IMAG
} real_imag_t;

/* own type to set the MPI thread level when MPI alone
 * or MPI and QMP are in use 
 * the default type is defined so that the type always exists, even
 * if the code is compiled without MPI support */

#ifdef QPHIX_QMP_COMMS
typedef enum tm_mpi_thread_level_t {
  TM_MPI_THREAD_SINGLE = QMP_THREAD_SINGLE,
  TM_MPI_THREAD_MULTIPLE = QMP_THREAD_MULTIPLE
} tm_mpi_thread_level_t;
#elif TM_USE_MPI
typedef enum tm_mpi_thread_level_t {
  TM_MPI_THREAD_SINGLE = MPI_THREAD_SERIALIZED,
  TM_MPI_THREAD_MULTIPLE = MPI_THREAD_MULTIPLE
} tm_mpi_thread_level_t;
#else
typedef enum tm_mpi_thread_level_t {
  TM_MPI_THREAD_SINGLE = 0,
  TM_MPI_THREAD_MULTIPLE
} tm_mpi_thread_level_t;
#endif

static inline void print_mpi_thread_level(const tm_mpi_thread_level_t thread_level)
{
  const char * const thread_level_string_single = "single";
  const char * const thread_level_string_multiple = "multiple";
  const char * selected_thread_level_string = thread_level_string_single;
  if( thread_level == TM_MPI_THREAD_MULTIPLE ){
    selected_thread_level_string = thread_level_string_multiple;
  }
  tm_debug_printf(0, 0, "MPI thread level set to '%s'\n", selected_thread_level_string);
}

/* types to track the state of a gauge field, could be added to a compound
 * gauge field type to provide meta-information or it can be kept
 * as a stand-alone object at the global scope to keep state */

typedef struct tm_GaugeHaloState_t {
  long long int gauge_id;
  int exchanged;
} tm_GaugeHaloState_t;

typedef struct tm_GaugeState_t {
  double gauge_id;
  int loaded;
  char name[TM_GAUGE_FIELD_NAME_LENGTH];
  tm_GaugeHaloState_t halo_state;
} tm_GaugeState_t;

static inline tm_GaugeState_t new_tm_GaugeState(const char * const name) {
  tm_GaugeState_t ret;
  ret.loaded = 0;
  ret.halo_state.exchanged = 0;
  ret.gauge_id = -TM_GAUGE_PROPAGATE_THRESHOLD;
  snprintf(ret.name, TM_GAUGE_FIELD_NAME_LENGTH, "%s", name);
  return(ret);
}

static inline void update_tm_gauge_id(tm_GaugeState_t * gauge_state, const double step){
  if(gauge_state->loaded){
    gauge_state->gauge_id += step;
    gauge_state->halo_state.exchanged = 0;
    tm_debug_printf(0, 4, 
                    "# [update_tm_gauge_id]: gauge id of %s stepped by %.12f to %.12f\n",
                    gauge_state->name, step, gauge_state->gauge_id);
  } else {
    gauge_state->loaded = 1;
    gauge_state->gauge_id = 0;
    gauge_state->halo_state.exchanged = 0;
    tm_debug_printf(0, 4, 
                    "# [update_tm_gauge_id]: gauge id of %s reset to 0\n", 
                    gauge_state->name);
  }
}

static inline void update_tm_gauge_exchange(tm_GaugeState_t * gauge_state){
  if(gauge_state->loaded){
    gauge_state->halo_state.gauge_id = gauge_state->gauge_id;
    gauge_state->halo_state.exchanged = 1;
    tm_debug_printf(0, 4, 
                    "# [update_tm_gauge_exchange]: %s tagged as freshly exchanged\n", 
                    gauge_state->name);
  }
}

static inline int check_tm_gauge_exchange(tm_GaugeState_t * gauge_state){
  return( gauge_state->halo_state.gauge_id == gauge_state->gauge_id && 
          gauge_state->halo_state.exchanged );
}

typedef struct tm_CloverState_t {
  long long int gauge_id;
  double c_sw;
  double kappa;
  int loaded;
} tm_CloverState_t;

static inline tm_CloverState_t new_tm_CloverState(void) {
  tm_CloverState_t ret;
  ret.loaded = 0;
  ret.gauge_id = -TM_GAUGE_PROPAGATE_THRESHOLD;
  return(ret);
}

typedef enum tm_CloverInverseType_t {
  // default inverse of the clover term, see sw_invert 
  TM_CLOVER_INVERSE_DEFAULT = 0, 
  // inverse of the clover term used for the ND operator, see sw_invert_nd
  TM_CLOVER_INVERSE_ND,
  // inverse of the ND clover term, multiplied by mubar, see sw_invert_mubar
  // this is used for convenience in the QPhiX interface
  TM_CLOVER_INVERSE_MUBAR,
  // inverse of the ND clover term, multiplied by epsbar, see sw_invert_epsbar
  // this is used for convenience in the QPhiX interface
  TM_CLOVER_INVERSE_EPSBAR
} tm_CloverInverseType_t;

typedef struct tm_CloverInverseState_t {
  tm_CloverInverseType_t inverse_type;
  long long int gauge_id;
  double c_sw;
  double kappa;
  double mu;
  double mubar;
  double epsbar;
  int loaded;
} tm_CloverInverseState_t;

static inline tm_CloverInverseState_t new_tm_CloverInverseState(void) {
  tm_CloverInverseState_t ret;
  ret.loaded = 0;
  ret.gauge_id = -TM_GAUGE_PROPAGATE_THRESHOLD;
  return(ret);
}

#endif // MISC_TYPES_H
