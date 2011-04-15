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


#include "spinor.ih"

int write_spinor(WRITER * writer, spinor ** const s, spinor ** const r, const int flavours, const int prec)
{
  DML_Checksum checksum;
  uint64_t bytes;
  int i = 0, status = 0;

  bytes = (n_uint64_t)LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * (n_uint64_t)(sizeof(spinor) * prec / 64);

  if(r == NULL) {
    for (i = 0; i < flavours; ++i) {
      //DEBUG following line
      write_header(writer, 1, 0, "scidac-binary-data", bytes);
      status  = write_binary_spinor_data_l(s[i], writer, &checksum, prec);
      write_checksum(writer, &checksum, NULL);
    }
  }
  else {
    for (i = 0; i < flavours; ++i) {
      write_header(writer, 1, 0, "scidac-binary-data", bytes);
      status = write_binary_spinor_data(s[i], r[i], writer, &checksum, prec);
      write_checksum(writer, &checksum, NULL);
    }
  }
  return status;
}
