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

int read_lemon_gauge_field_parallel(char *filename)
{
  MPI_File ifs;
  int status;
  char * header_type;
  LemonReader * lemonreader;
  DML_Checksum checksum_read;
  DML_Checksum checksum_calc;
  int DML_read_flag = 0;
  int gauge_read_flag = 0;

  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &ifs);
  if (status)
  {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonreader = lemonCreateReader(&ifs, g_cart_grid);
  if( lemonreader == (LemonReader *)NULL )
  {
    fprintf(stderr, "Unable to open LemonReader\n");
    MPI_File_close(&ifs);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  while ((status = lemonReaderNextRecord(lemonreader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "lemonReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = lemonReaderType(lemonreader);

    if (strcmp("ildg-binary-data", header_type) == 0)
    {
      read_binary_gauge_data_parallel(lemonreader, &checksum_calc);
      gauge_read_flag = 1;
    }
    else if (strcmp("scidac-checksum", header_type) == 0)
    {
      if (g_debug_level > 1)
      {
        read_checksum_parallel(lemonreader, &checksum_read);
        DML_read_flag = 1;
      }
    }
    lemonReaderCloseRecord(lemonreader);
  }

  if (!gauge_read_flag)
  {
    if (g_cart_id == 0)
    {
      fprintf(stderr, "Did not find gauge record in %s.\n", filename);
      fprintf(stderr, "Panic! Aborting...\n");
      fflush(stderr);
    }
    lemonDestroyReader(lemonreader);
    MPI_File_close(&ifs);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(501);
  }

  if (g_debug_level > 1 && g_cart_id == 0)
  {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# No DML checksum record found in %s.\n", filename);
    fflush(stderr);
  }

  lemonDestroyReader(lemonreader);
  MPI_File_close(&ifs);
  return 0;
}
