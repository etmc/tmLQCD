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

int write_message(WRITER * writer, char const *buffer, MPI_Offset bytes)
{
  int status;
  MPI_Offset bytesWritten = bytes;

#ifndef HAVE_LIBLEMON
  if(g_cart_id == 0){
#endif /* ! HAVE_LIBLEMON */
    if (buffer == (char*)NULL)
      return(0);

    status = WriteRecordData((void*)buffer, &bytes, writer);
    if (status != LIME_SUCCESS || bytes != bytesWritten)
      kill_with_error(writer->fp, g_cart_id, "I/O error on writing message. Aborting...\n");
#ifndef HAVE_LIBLEMON
  }
#endif /* ! HAVE_LIBLEMON */
  return(0);
}
