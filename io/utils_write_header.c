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

void write_header(LimeWriter * limewriter, int MB, int ME, char *type, uint64_t bytes)
{
  int status;
  LimeRecordHeader *limeheader;

  limeheader = limeCreateHeader(MB, ME, type, bytes);
  status = limeWriteRecordHeader(limeheader, limewriter);
  limeDestroyHeader(limeheader);

  if (status != LIME_SUCCESS)
    kill_with_error(limewriter->fp, g_cart_id, "LEMON header writing error. Aborting\n");
  return;
}
