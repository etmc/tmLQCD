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
#define TM_GAUGE_PROPAGATE_THRESHOLD 10.0
#define TM_GAUGE_PROPAGATE_MIN 0.01

/* enumeration type for the identity of the program
 * which is being executed
 * this is useful to unify various utility functions which
 * otherwise lead to a lot of code duplication */
typedef enum tm_ProgramId_t {
  TM_PROGRAM_HMC_TM = 0,
  TM_PROGRAM_INVERT,
  TM_PROGRAM_OFFLINE_MEASUREMENT,
  TM_PROGRAM_BENCHMARK,
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

typedef enum backup_restore_t {
  TM_BACKUP_GLOBALS = 0,
  TM_RESTORE_GLOBALS
} backup_restore_t;

typedef enum real_imag_t {
  TM_REAL = 0,
  TM_IMAG
} real_imag_t;

#endif // MISC_TYPES_H
