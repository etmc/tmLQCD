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
  int status = 0, getpos = 0, bytes = 0, prec = 0, prop_type, position = position_;
  char *header_type = NULL;
  READER *reader = NULL;
  DML_Checksum checksum;

  construct_reader(&reader, filename);
  /* determine the propagator type */

  prop_type = parse_propagator_type(reader);
  if ( prop_type == -1 ) {
    kill_with_error(reader->fp, g_cart_id, "Did not find propagator or source type. Aborting...\n");
/*     if(g_proc_id == 0) fprintf(stderr, "hack here, please prepare sources with correct headers!\n"); */
/*     prop_type = 10; */
  }

  /* strictly speeking the following depends on whether we read a source or a propagator */
  if(prop_type == 1) position = 2*position_ + 1;
  /* anything else needs implementation! */
  else if(prop_type == 2 || prop_type == 3) 
    kill_with_error(reader->fp, g_cart_id, "Propagator type not yet implemented. Aborting read!\n");
  else if(prop_type == 11 || prop_type == 12 || prop_type == 13) 
    kill_with_error(reader->fp, g_cart_id, "Source type not yet implemented. Aborting read!\n");
  else if(prop_type == -1)
    kill_with_error(reader->fp, g_cart_id, "No propagator or source type record in file. Aborting read!\n");

  /* Find the desired propagator (could be more than one in a file) */
  while ((status = ReaderNextRecord(reader)) != LIME_EOF) {
    if (status != LIME_SUCCESS) {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = ReaderType(reader);
    if (strcmp("scidac-binary-data", header_type) == 0) {
      if (getpos == position)
        break;
      else
        ++getpos;
    }
  }

  if (status == LIME_EOF)
    kill_with_error(reader->fp, g_cart_id, "No scidac-binary-data record found in file.\n");

  bytes = ReaderBytes(reader);
  if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor))
    prec = 64;
  else
    if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor) / 2)
      prec = 32;
    else
      kill_with_error(reader->fp, g_cart_id, "Wrong length in eospinor. Aborting read!\n");

  if (g_cart_id == 0 && g_debug_level > 2)
    printf("# %d bit precision read.\n", prec);

  read_binary_spinor_data(s, r, reader, &checksum);

  if (g_cart_id == 0 && g_debug_level > 1)
    printf("# checksum for DiracFermion field in file %s position %d is %#x %#x\n", filename, position, checksum.suma, checksum.sumb);

  destruct_reader(reader);

  return 0;
}
