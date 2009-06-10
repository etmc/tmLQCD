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

void write_header_parallel(LemonWriter * lemonwriter, int MB, int ME, char *type, uint64_t bytes)
{
  int status;
  LemonRecordHeader *lemonheader;

  lemonheader = lemonCreateHeader(MB, ME, type, bytes);
  status = lemonWriteRecordHeader(lemonheader, lemonwriter);
  lemonDestroyHeader(lemonheader);

  if (status != LEMON_SUCCESS)
    kill_with_error(lemonwriter->fh, lemonwriter->my_rank, "LEMON header writing error. Aborting\n");
}
