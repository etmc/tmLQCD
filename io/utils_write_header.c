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

void write_header(WRITER * writer, int MB, int ME, char const *type, uint64_t bytes)
{
  int status;
  RECORD_HEADER *header;

#ifndef HAVE_LIBLEMON
  if(g_cart_id == 0) {
#endif /* ! HAVE_LIBLEMON */
    /* Nasty (but probably harmless) hack to get rid of const qualifier - the original c-lime was sloppy here. */
    header = CreateHeader(MB, ME, (char*)type, bytes);
    status = WriteRecordHeader(header, writer);
    DestroyHeader(header);

    if (status != LIME_SUCCESS) {
      kill_with_error(writer->fp, g_cart_id, "Header writing error. Aborting\n");
    }
#ifndef HAVE_LIBLEMON
  }
#endif /* ! HAVE_LIBLEMON */
  return;
}
