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

#include "utils.ih"

int write_checksum_parallel(LemonWriter * lemonwriter, DML_Checksum * checksum)
{
  LemonRecordHeader * lemonheader = NULL;
  int status = 0;
  int MB_flag = 0, ME_flag = 1;
  char *message;
  uint64_t bytes;

  message = (char*)malloc(512);
  if (message == (char*)NULL )
  {
    fprintf(stderr, "Memory error in write_checksum_parallel. Aborting\n");
    MPI_File_close(lemonwriter->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                   "<scidacChecksum>\n"
                   "  <version>1.0</version>\n"
                   "  <suma>%x</suma>\n"
                   "  <sumb>%x</sumb>\n"
                   "</scidacChecksum>", checksum->suma, checksum->sumb);
  bytes = strlen( message );
  lemonheader = lemonCreateHeader(MB_flag, ME_flag, "scidac-checksum", bytes);
  status = lemonWriteRecordHeader( lemonheader, lemonwriter);
  if(status < 0 ) {
    fprintf(stderr, "LEMON write header error %d\n", status);
    MPI_File_close(lemonwriter->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader( lemonheader );
  lemonWriteRecordData(message, &bytes, lemonwriter);
  lemonWriterCloseRecord(lemonwriter);
  free(message);
  return(0);
}
