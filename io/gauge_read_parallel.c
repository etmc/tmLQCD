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

void read_lemon_gauge_field_parallel(char *filename, char **scidac_checksum,
			 char **xlf_info, char **ildg_data_lfn)
{
  MPI_File ifs;
  int status;
  char *header_type = NULL;
  LemonReader *lemonreader = NULL;
  DML_Checksum checksum_calc;
  DML_Checksum checksum_read;
  int DML_read_flag = 0;
  int gauge_read_flag = 0;
  char *checksum_string = NULL;

  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &ifs);
  if (status != MPI_SUCCESS)
    kill_with_error(&ifs, g_cart_id, "Could not open file. Aborting...\n");

  lemonreader = lemonCreateReader(&ifs, g_cart_grid);
  if (lemonreader == (LemonReader *)NULL)
    kill_with_error(&ifs, g_cart_id, "Could not create lemon reader. Aborting...\n");

  while ((status = lemonReaderNextRecord(lemonreader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "lemonReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = lemonReaderType(lemonreader);
    fprintf(stdout, "Read header: %s\n", header_type); fflush(stdout);

    if (strcmp("ildg-binary-data", header_type) == 0)
    {
      read_binary_gauge_data_parallel(lemonreader, &checksum_calc);
      gauge_read_flag = 1;
    }
    else if (strcmp("scidac-checksum", header_type) == 0)
    {
      checksum_string = (char*)malloc(1024);
      read_message_parallel(lemonreader, &checksum_string);
      checksum_string[1023] = '\0';
      DML_read_flag = parse_checksum_xml(checksum_string, &checksum_read);
      if (scidac_checksum != (char**)NULL)
      {
        if ((*scidac_checksum) != NULL)
          free(*scidac_checksum);
        *scidac_checksum = (char*)malloc(strlen(checksum_string) + 1);
        strcpy(*scidac_checksum, checksum_string);
      }
      free(checksum_string);
    }
    else if (strcmp("xlf-info", header_type) == 0)
    {
      read_message_parallel(lemonreader, xlf_info);
    }
    else if (strcmp("ildg-data-lfn", header_type) == 0)
    {
      read_message_parallel(lemonreader, ildg_data_lfn);
    }
    lemonReaderCloseRecord(lemonreader);
  }
  
  if (!gauge_read_flag)
    kill_with_error(&ifs, g_cart_id, "Did not find gauge record. Aborting...\n");

  if (g_debug_level > 1 && g_cart_id == 0)
  {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# Scidac checksum record not found or not read.\n");
    fflush(stdout);
  }

  fprintf(stdout, "Ready to destroy reader.\n"); fflush(stdout);
  lemonDestroyReader(lemonreader);
  MPI_File_close(&ifs);
}
