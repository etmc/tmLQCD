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

#endif // MISC_TYPES_H
