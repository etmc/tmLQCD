/***********************************************************************
* Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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

#include "spinor.ih"
#include "default_input_values.h"

paramsPropInfo PropInfo = {_default_propagator_splitted, _default_source_format_flag, _default_prop_precision_flag, NULL};
paramsSourceInfo SourceInfo = {0, _default_propagator_splitted, _default_source_format_flag, _default_prop_precision_flag, 0, 0, 0, 0, 0, 0, 0, 1, NULL};

int read_spinor(spinor * const s, spinor * const r, char * filename, const int position_) {
  int status = 0, getpos = 0, bytes = 0, prec = 0, prop_type, position = position_, rstat=0;
  char *header_type = NULL;
  READER *reader = NULL;
  DML_Checksum checksum;
  construct_reader(&reader, filename);
  /* determine the propagator type */
  prop_type = parse_propagator_type(reader);

  switch (prop_type) {
  case 1:
    /* strictly speaking the following depends on whether we read a source or a propagator */
    position = 2 * position_ +1;
    break;
  case 2:
  case 3:
    return(-2);
  case 11:
  case 12:
  case 13:
    return(-3);
  case -1:
  case 4:
    prop_type = 0;
  }

  /* seek back to beginning of file*/
  destruct_reader(reader);
  construct_reader(&reader, filename);

  /* Find the desired propagator (could be more than one in a file) */
  while ((status = ReaderNextRecord(reader)) != LIME_EOF) {
    if (status != LIME_SUCCESS) {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = ReaderType(reader);
    if (strcmp("scidac-binary-data", header_type) == 0) {
      if (getpos == position) {
        break;
      }
      ++getpos;
    }
  }

  if (status == LIME_EOF) {
    fprintf(stderr, "Unable to find requested LIME record scidac-binary-data in file %s.\nEnd of file reached before record was found.\n", filename);
    return(-5);
  }

  bytes = ReaderBytes(reader);

  if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor)) {
    prec = 64;
  }
  else {
    if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor) / 2) {
      prec = 32;
    }
    else {
      fprintf(stderr, "Length of scidac-binary-data record in %s does not match input parameters.\n", filename);
      fprintf(stderr, "Found %d bytes.\n", bytes);
      return(-6);
    }
  }

  if (g_cart_id == 0 && g_debug_level >= 0) {
    printf("# %s precision read (%d bits).\n", (prec == 64 ? "Double" : "Single") ,prec);
  }

  if(r == NULL) {
    if( (rstat = read_binary_spinor_data_l(s, reader, &checksum)) != 0) {
      fprintf(stderr, "read_binary_spinor_data_l failed with return value %d", rstat);
      return(-7);
    }
  }
  else {
    if( (rstat = read_binary_spinor_data(s, r, reader, &checksum)) != 0) {
      fprintf(stderr, "read_binary_spinor_data failed with return value %d", rstat);
      return(-7);
    }
  }

  if (g_cart_id == 0 && g_debug_level >= 0) {
    printf("# Scidac checksums for DiracFermion field %s position %d:\n", filename, position);
    printf("#   Calculated            : A = %#x B = %#x.\n", checksum.suma, checksum.sumb);
    printf("# No Scidac checksum was read from headers, unable to check integrity of file.\n");
  }

  destruct_reader(reader);

  return(0);
}
