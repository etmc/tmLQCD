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

int read_message(READER * reader, char **buffer)
{
  int status;
  n_uint64_t bytes, read_bytes;

#ifndef HAVE_LIBLEMON
  if (g_cart_id == 0){
#endif /* HAVE_LIBLEMON */
  if (buffer == (char**)NULL)
    return(-1);
    
  if ((*buffer) != (char*)NULL)
    free(*buffer);

  bytes = ReaderBytes(reader);

  *buffer = (char*)malloc(bytes + 1);

  if (*buffer  == (char*)NULL)
  {
    fprintf(stderr, "Couldn't malloc data buffer in read_message.\n");
    return(-1);
  }

  status = ReaderReadData(*buffer, &bytesRead, reader);

  if (status != LEMON_SUCCESS || bytes != bytesRead)
    kill_with_error(reader->fp, g_cart_id, "Error in reading message.\n");

  buffer[bytes] = '\0'; /* Force termination for safety */
#ifndef HAVE_LIBLEMON
  }
#endif /* HAVE_LIBLEMON */
  return(0);
}
