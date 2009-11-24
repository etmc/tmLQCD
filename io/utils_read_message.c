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

int read_message(LimeReader * limereader, char **buffer) {
  
  int status;
  n_uint64_t bytes, read_bytes;

  if (buffer == (char**)NULL)
    return(-1);
  
  if ((*buffer) != (char*)NULL)
    free(*buffer);


  bytes = limeReaderBytes(limereader);

  if((*buffer = (char*)malloc(bytes + 1)) == (char*)NULL) {
    fprintf(stderr, "Couldn't malloc data buf in read_message\n");
    return(-1);
  }
  read_bytes = bytes;
  status = limeReaderReadData((void *)*buffer, &read_bytes, limereader);

  if( status < 0 ) {
    if( status != LIME_EOR ) {
      fprintf(stderr, "LIME read error occurred: status= %d  %llu bytes wanted, %llu read\n",
	      status, (unsigned long long)bytes,
	      (unsigned long long)read_bytes);
      free(*buffer);
      *buffer = NULL;
      return(-1);
    }
  }
  buffer[bytes]='\0';
  return(0);
}
