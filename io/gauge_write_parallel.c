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

int write_lemon_gauge_field_parallel(char * filename, const double plaq, const int counter, const int prec) {
  MPI_File ofs;
  LemonWriter * lemonwriter = NULL;
  LemonRecordHeader * lemonheader = NULL;

  int status=0;
  uint64_t bytes;
  DML_Checksum checksum;

  MPI_File_open(g_cart_grid, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ofs);

  lemonwriter = lemonCreateWriter( &ofs, g_cart_grid );
  if(lemonwriter == (LemonWriter*)NULL) {
    fprintf(stderr, "LEMON error in file %s for writing!\nPanic! Aborting...\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  write_xlf_info_parallel(lemonwriter, plaq, counter);

  bytes = ((uint64_t)LX*g_nproc_x)*((uint64_t)LY*g_nproc_y)*((uint64_t)LZ*g_nproc_z)*((uint64_t)T*g_nproc_t)*((uint64_t)4*sizeof(su3));
  if(prec == 32) bytes = bytes/((uint64_t)2);

  write_ildg_format_parallel(lemonwriter, prec);

  lemonheader = lemonCreateHeader(0, 0, "ildg-binary-data", bytes);
  status = lemonWriteRecordHeader(lemonheader, lemonwriter);
  if(status < 0 ) {
    fprintf(stderr, "LEMON write header error %d\n", status); fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader(lemonheader);

  write_binary_gauge_data_parallel(lemonwriter, prec, &checksum);

  if(g_cart_id == 0) {
    fprintf(stderr, "Checksum A: %#x \nChecksum B: %#x\n", checksum.suma, checksum.sumb);
    fflush(stderr);
  }
  write_checksum_parallel(lemonwriter, &checksum);
  MPI_Barrier(MPI_COMM_WORLD);
  lemonDestroyWriter( lemonwriter );
  MPI_File_close(&ofs);
  return(0);
}
