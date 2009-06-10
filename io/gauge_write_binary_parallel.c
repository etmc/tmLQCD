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

#include "gauge.ih"

void write_binary_gauge_data_parallel(LemonWriter * lemonwriter, const int prec, DML_Checksum * ans)
{
  int x, xG, y, yG, z, zG, t, tG;
  su3 tmp3[4];
  int globaldims[] = {L, L, L, T_global};
  unsigned long bufoffset;
  char *filebuffer;
  uint64_t bytes;
  DML_SiteRank rank;
  DML_checksum_init(ans);
  bytes = (uint64_t)sizeof(su3) * (prec == 32 ? 2 : 4);
  bufoffset = 0;
  filebuffer = (char*)malloc(bytes * VOLUME);
  double tick = 0, tock = 0;
  char measure[64];

  tG = g_proc_coords[0]*T;
  zG = g_proc_coords[3]*LZ;
  yG = g_proc_coords[2]*LY;
  xG = g_proc_coords[1]*LX;
  for(t = 0; t < T; t++) {
      for(z = 0; z < LZ; z++) {
        for(y = 0; y < LY; y++) {
          for(x = 0; x < LX; x++) {
            rank = (DML_SiteRank) ((((tG + t)*L + zG + z)*L + yG + y)*L + xG + x);
            memcpy(&tmp3[0], &g_gauge_field[ g_ipt[t][x][y][z] ][1], sizeof(su3));
            memcpy(&tmp3[1], &g_gauge_field[ g_ipt[t][x][y][z] ][2], sizeof(su3));
            memcpy(&tmp3[2], &g_gauge_field[ g_ipt[t][x][y][z] ][3], sizeof(su3));
            memcpy(&tmp3[3], &g_gauge_field[ g_ipt[t][x][y][z] ][0], sizeof(su3));
  #ifndef WORDS_BIGENDIAN
            if(prec == 32)
              byte_swap_assign_double2single(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            else
              byte_swap_assign(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
  #else
            if(prec == 32)
              double2single(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            else
              memcpy(filebuffer + bufoffset, tmp3, 4*sizeof(su3));
  #endif
            DML_checksum_accum(ans, rank, (char*) filebuffer + bufoffset, bytes);

            bufoffset += bytes;
          }
        }
      }
    }


  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  lemonWriteLatticeParallel(lemonwriter, filebuffer, bytes, globaldims);

  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0)
    {
      engineering(measure, L * L * L * T_global * bytes, "b");
      fprintf(stdout, "Time spent writing %s ", measure);
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (tock - tick), "b/s");
      fprintf(stdout, "Writing speed: %s", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }

  lemonWriterCloseRecord(lemonwriter);

  DML_global_xor(&ans->suma);
  DML_global_xor(&ans->sumb);

  free(filebuffer);
}
