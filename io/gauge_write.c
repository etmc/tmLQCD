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

void write_gauge_field(char * filename, const int prec, paramsXlfInfo const *xlfInfo)
{
  WRITER * writer = NULL;
  LIME_FILE *ofs = NULL;
  uint64_t bytes;

  DML_Checksum     checksum;
  paramsIldgFormat *ildg;

#ifdef HAVE_LIBLEMON
  MPI_File ofs_object;

  ofs = &ofs_object;
  MPI_File_open(g_cart_grid, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ofs);
#else /* HAVE_LIBLEMON */
  if(g_cart_id == 0)
    ofs = fopen(filename, "w");
#endif /* HAVE_LIBLEMON */

  bytes = (uint64_t)L * L * L * T_global * sizeof(su3) * prec / 16;

  writer = CreateWriter(ofs, g_cart_grid);
  if (writer == (WRITER*)NULL)
#ifdef HAVE_LIBLEMON
    kill_with_error(writer->fh, writer->my_rank, "Cannot construct writer. Aborting...\n");
#else /* HAVE_LIBLEMON */
    kill_with_error(writer->fp, writer->my_rank, "Cannot construct writer. Aborting...\n");
#endif /* HAVE_LIBLEMON */

  write_xlf_info(writer, xlfInfo);

  ildg = construct_paramsIldgFormat(prec);
  write_ildg_format(writer, ildg);
  free(ildg);

  write_header(writer, 0, 0, "ildg-binary-data", bytes);
  write_binary_gauge_data(writer, prec, &checksum);
  write_checksum(writer, &checksum, NULL);

  if (g_cart_id == 0)
  {
    fprintf(stdout, "# Checksum A: %#x \nChecksum B: %#x\n", checksum.suma, checksum.sumb);
    fflush(stdout);
#ifndef HAVE_LIBLEMON
    limeDestroyWriter(writer);
    fclose(ofs);
#endif /* ! HAVE_LIBLEMON */
  }
#ifdef HAVE_LIBLEMON
  MPI_Barrier(MPI_COMM_WORLD);
  DestroyWriter(writer);
  MPI_File_close(&ofs);
#endif /* HAVE_LIBLEMON */

  return;
}
