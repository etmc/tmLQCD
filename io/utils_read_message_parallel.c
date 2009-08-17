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

void read_message_parallel(LemonReader * lemonreader, char **buffer)
{
  int status;
  uint64_t bytes = lemonReaderBytes(lemonreader);
  uint64_t bytesRead = bytes;

  if (buffer == (char**)NULL)
    return;
  
  if ((*buffer) != (char*)NULL)
    free(*buffer);

  *buffer = (char*)malloc(bytes + 1);
  buffer[bytes] = '\0'; /* Force termination for safety */
  status = lemonReaderReadData(*buffer, &bytesRead, lemonreader);

  if (status != LEMON_SUCCESS || bytes != bytesRead)
    if (lemonreader->my_rank == 0)
    {
      fprintf(stderr, "Error in reading message.\n");
      MPI_File_close(lemonreader->fh);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(500);
    }
}
