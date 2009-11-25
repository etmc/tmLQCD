/***********************************************************************
* Copyright (C) 2009 Siebren Recker, Albert Deuzeman, Carsten Urbach
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

#include "spinor.ih"

void read_spinor_parallel(spinor * const s, spinor * const r,
                          char * filename, const int position)
{
  MPI_File ifs;
  int status = 0, getpos = 0;
  uint64_t bytes;
  char *header_type;
  LemonReader *lemonreader;
  int prec = 32;
  DML_Checksum checksum;

  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &ifs);
  if (status)
    kill_with_error(&ifs, g_cart_id, "Unable to open file.\n");

  lemonreader = lemonCreateReader(&ifs, g_cart_grid);
  if (lemonreader == (LemonReader *)NULL)
    kill_with_error(&ifs, g_cart_id, "Unable to create lemon reader.\n");

  /* Find the desired propagator (could be more than one in a file) */
  while ((status = lemonReaderNextRecord(lemonreader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "lemonReaderNextRecord returned status %d.\n", status);
      break;
    }
    header_type = lemonReaderType(lemonreader);

    if (strcmp("scidac-binary-data", header_type) == 0)
    {
      if (getpos == position)
        break;
      else
        ++getpos;
    }
  }

  if (status == LEMON_EOF)
    kill_with_error(&ifs, g_cart_id, "No scidac-binary-data record found in file.\n");

  bytes = lemonReaderBytes(lemonreader);
  if ((int)bytes == LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor))
    prec = 64;
  else
    if ((int)bytes == LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor) / 2)
      prec = 32;
    else
      kill_with_error(&ifs, g_cart_id, "Wrong length in eospinor. Aborting read!\n");
  if (g_cart_id == 0 && g_debug_level > 2)
    printf("# %d bit precision read\n", prec);

  read_binary_spinor_data_parallel(s, r, lemonreader, &checksum);

  if (g_cart_id == 0 && g_debug_level > 1)
    printf("# checksum for DiracFermion field in file %s position %d is %#x %#x\n",
           filename, position, checksum.suma, checksum.sumb);

  lemonDestroyReader(lemonreader);
  MPI_File_close(&ifs);
}
