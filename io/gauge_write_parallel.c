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

void write_lemon_gauge_field_parallel(char * filename, int prec, paramsXlfInfo const *xlfInfo)
{
  MPI_File ofs;
  LemonWriter * writer = NULL;

  uint64_t bytes;

  DML_Checksum     checksum;
  paramsIldgFormat *ildg;

  bytes = (uint64_t)L * L * L * T * sizeof(su3) * prec / 16;

  MPI_File_open(g_cart_grid, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ofs);

  writer = lemonCreateWriter(&ofs, g_cart_grid);
  if (writer == (LemonWriter*)NULL)
    kill_with_error(writer->fh, writer->my_rank, "Cannot construct LemonWriter. Aborting...\n");

  write_xlf_info_parallel(writer, xlfInfo);

  ildg = construct_paramsIldgFormat(prec);
  write_ildg_format_parallel(writer, ildg);
  free(ildg);

  write_header_parallel(writer, 0, 0, "ildg-binary-data", bytes);
  write_binary_gauge_data_parallel(writer, prec, &checksum);
  write_checksum_parallel(writer, &checksum);

  if (g_cart_id == 0)
  {
    fprintf(stdout, "Checksum A: %#x \nChecksum B: %#x\n", checksum.suma, checksum.sumb);
    fflush(stdout);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  lemonDestroyWriter(writer);
  MPI_File_close(&ofs);
}
