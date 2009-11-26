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

#include "gauge.ih"

int read_gauge_field(char * filename, DML_Checksum *scidac_checksum,
                     char **xlf_info, char **ildg_data_lfn) {
  int status = 0;
  char *header_type = NULL;
  READER *reader = NULL;

  DML_Checksum checksum_read;
  DML_Checksum checksum_calc;
  int DML_read_flag = 0;
  int gauge_read_flag = 0;
  char *checksum_string = NULL;
  LIME_FILE *ifs = NULL;

#ifdef HAVE_LIBLEMON
  MPI_File ifs_object;

  ifs = &ifs_object;
  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, ifs);
  if (status != MPI_SUCCESS)
    kill_with_error(ifs, g_cart_id, "Could not open file. Aborting...\n");
#else /* HAVE_LIBLEMON */
  ifs = fopen(filename, "r");
  if (ifs == (FILE *)NULL)
    kill_with_error(ifs, g_cart_id, "Could not open file. Aborting...\n");
#endif /* HAVE_LIBLEMON */

  reader = CreateReader(&ifs, g_cart_grid);
  if (reader == (READER *)NULL)
    kill_with_error(ifs, g_cart_id, "Could not create reader. Aborting...\n");

  while ((status = ReaderNextRecord(reader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = ReaderType(reader);
    if(g_cart_id == 0 && g_debug_level > 1) {
      fprintf(stderr, "found header %s, will now read the message\n", header_type);
    }

    if (strcmp("ildg-binary-data", header_type) == 0)
    {
      read_binary_gauge_data(reader, &checksum_calc);
      gauge_read_flag = 1;
    }
    else if (strcmp("scidac-checksum", header_type) == 0)
    {
      read_message(reader, &checksum_string);
      DML_read_flag = parse_checksum_xml(checksum_string, &checksum_read);
      if (DML_read_flag && scidac_checksum != (DML_Checksum*)NULL)
        *scidac_checksum = checksum_read;
      free(checksum_string);
    }
    else if (strcmp("xlf-info", header_type) == 0)
    {
      read_message(reader, xlf_info);
    }
    else if (strcmp("ildg-data-lfn", header_type) == 0)
    {
      read_message(reader, ildg_data_lfn);
    }
    ReaderCloseRecord(reader);
  }

  if (!gauge_read_flag)
    kill_with_error(&ifs, g_cart_id, "Did not find gauge record. Aborting...\n");


  if (g_debug_level > 0 && g_cart_id == 0)
  {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# Scidac checksum record not found or malformed.\n");
    fflush(stdout);
  }

  DestroyReader(reader);

#ifdef HAVE_LIBLEMON
  MPI_File_close(ifs);
#else /* HAVE_LIBLEMON */
  fclose(ifs);
#endif /* HAVE_LIBLEMON */

  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;

  return 0;
}
