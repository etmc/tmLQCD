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

void read_spinor(spinor * const s, spinor * const r,
		 char * filename, const int position)
{
#ifdef HAVE_LIBLEMON
  read_spinor_parallel(s, r, filename, position);
#else
  FILE * ifs = NULL;
  int status = 0, getpos = 0;
  char *header_type;
  LimeReader *limereader;
  DML_Checksum checksum;

  ifs = fopen(filename, "r");
  if (ifs == NULL)
    kill_with_error(ifs, g_cart_id, "Unable to open file.\n");

  limereader = limeCreateReader(ifs);
  if (limereader == (LimeReader *)NULL)
    kill_with_error(ifs, g_cart_id, "Unable to create lime reader.\n");

  /* Find the desired propagator (could be more than one in a file) */
  while ((status = limeReaderNextRecord(limereader)) != LIME_EOF) {

    if (status != LIME_SUCCESS) {
      fprintf(stderr, "limeReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = limeReaderType(limereader);
    if (strcmp("scidac-binary-data", header_type) == 0) {
      if (getpos == position)
        break;
      else
        ++getpos;
    }
  }

  if (status == LIME_EOF)
    kill_with_error(ifs, g_cart_id, "No scidac-binary-data record found in file.\n");

  read_binary_spinor_data(s, r, limereader, &checksum);

  if (g_cart_id == 0 && g_debug_level > 1)
    printf("# checksum for DiracFermion field in file %s position %d is %#x %#x\n",
           filename, position, checksum.suma, checksum.sumb);

  limeDestroyReader(limereader);
  fclose(ifs);
#endif
  return;
}
