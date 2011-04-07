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

/* FIXME I will first fix this function by using referral.
         Probably should be done better in the future. AD. */

#ifdef HAVE_LIBLEMON
int write_binary_gauge_data(LemonWriter * lemonwriter, const int prec, DML_Checksum * checksum)
{
  int x, xG, y, yG, z, zG, t, tG, status = 0;
  su3 tmp3[4];
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  int scidacMapping[] = {0, 3, 2, 1};
  unsigned long bufoffset;
  char * filebuffer = NULL;
  uint64_t bytes;
  double tick = 0, tock = 0;
  char measure[64];
  DML_SiteRank rank;
  DML_checksum_init(checksum);

  bytes = (uint64_t)sizeof(su3) * (prec == 32 ? 2 : 4);
  bufoffset = 0;
  if((void*)(filebuffer = (char*)malloc(bytes * VOLUME)) == NULL) {
    fprintf (stderr, "malloc errno in write_binary_gauge_data_parallel: %d\n",errno);
    fflush(stderr);
    errno = 0;
    return 1;
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

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
          if(prec == 32)
            be_to_cpu_assign_double2single(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
          else
            be_to_cpu_assign(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
          DML_checksum_accum(checksum, rank, (char*) filebuffer + bufoffset, bytes);
          bufoffset += bytes;
        }
      }
    }
  }

  status = lemonWriteLatticeParallelMapped(lemonwriter, filebuffer, bytes, latticeSize, scidacMapping);

  if (status != LEMON_SUCCESS)
  {
    free(filebuffer);
    fprintf(stderr, "LEMON write error occurred with status = %d, while writing in gauge_write_binary.c!\n", status);
    return(-2);
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes, "b");
      fprintf(stdout, "# Time spent writing %s ", measure);
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock - tick), "b/s");
      fprintf(stdout, "# Writing speed: %s", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }

  lemonWriterCloseRecord(lemonwriter);

  free(filebuffer);

  status = DML_global_xor(&checksum->suma);
  if (status != MPI_SUCCESS) {
    fprintf(stderr, "DML Checksum accumulation error occurred with status = %d, while writing in gauge_write_binary.c!\n", status);
    return(-2);
  }
  status = DML_global_xor(&checksum->sumb);
  if (status != MPI_SUCCESS) {
    fprintf(stderr, "DML Checksum accumulation error occurred with status = %d, while writing in gauge_write_binary.c!\n", status);
    return(-2);
  }

  return 0;
}

#else /* HAVE_LIBLEMON */

int write_binary_gauge_data(LimeWriter * limewriter, const int prec, DML_Checksum * checksum)
{
  int x, X, y, Y, z, Z, tt, t0, tag=0, id=0, status=0;
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  su3 tmp[4];
  su3 tmp3[4];
  float tmp2[72];
  int coords[4];
  n_uint64_t bytes;
  DML_SiteRank rank;
#ifdef MPI
  double tick = 0, tock = 0;
  char measure[64];
  MPI_Status mpi_status;
#endif

  DML_checksum_init(checksum);

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
#endif
  if(prec == 32) bytes = (n_uint64_t)2*sizeof(su3);
  else bytes = (n_uint64_t)4*sizeof(su3);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    tt = t0 - g_proc_coords[0]*T;
    coords[0] = t0 / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
        tag = 0;
        Y = y - g_proc_coords[2]*LY;
        coords[2] = y / LY;
        for(x = 0; x < LX*g_nproc_x; x++) {
          X = x - g_proc_coords[1]*LX;
          coords[1] = x / LX;
#ifdef MPI
          MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
          if(g_cart_id == 0) {
            /* Rank should be computed by proc 0 only */
            rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);
            if(g_cart_id == id) {
              memcpy(&tmp3[0], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], sizeof(su3));
              memcpy(&tmp3[1], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], sizeof(su3));
              memcpy(&tmp3[2], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], sizeof(su3));
              memcpy(&tmp3[3], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], sizeof(su3));

              if(prec == 32) {
                be_to_cpu_assign_double2single(tmp2, tmp3, 4*sizeof(su3)/8);
                DML_checksum_accum(checksum, rank, (char*) tmp2, 4*sizeof(su3)/2);
                status = limeWriteRecordData((void*)&tmp2, &bytes, limewriter);
              }
              else {
                be_to_cpu_assign(tmp, tmp3, 4*sizeof(su3)/8);
                DML_checksum_accum(checksum, rank, (char*) tmp, 4*sizeof(su3));
                status = limeWriteRecordData((void*)&tmp, &bytes, limewriter);
              }
            }
#ifdef MPI
            else {
              if(prec == 32) {
                MPI_Recv(tmp2, 4*sizeof(su3)/8, MPI_FLOAT, id, tag, g_cart_grid, &mpi_status);
                DML_checksum_accum(checksum, rank, (char*) tmp2, 4*sizeof(su3)/2);
                status = limeWriteRecordData((void*)&tmp2, &bytes, limewriter);
              }
              else {
                MPI_Recv(tmp, 4*sizeof(su3)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpi_status);
                DML_checksum_accum(checksum, rank, (char*) tmp, 4*sizeof(su3));
                status = limeWriteRecordData((void*)&tmp, &bytes, limewriter);
              }
            }
#endif
            if(status < 0 ) {
              fprintf(stderr, "LIME write error occurred with status = %d, while writing in gauge_write_binary.c!\n", status);
              fprintf(stderr, "x %d, y %d, z %d, t %d (%d,%d,%d,%d)\n",x,y,z,tt,X,Y,Z,tt);
              fprintf(stderr, "id = %d, bytes = %lu, size = %d\n", g_cart_id, bytes,  (int)(4*sizeof(su3)/8));
#ifdef MPI
              MPI_Abort(MPI_COMM_WORLD, 1);
              MPI_Finalize();
#endif
              exit(500);
            }
          }
#ifdef MPI
          else {
            if(g_cart_id == id){
              memcpy(&tmp3[0], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], sizeof(su3));
              memcpy(&tmp3[1], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], sizeof(su3));
              memcpy(&tmp3[2], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], sizeof(su3));
              memcpy(&tmp3[3], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], sizeof(su3));
              if(prec == 32) {
                be_to_cpu_assign_double2single(tmp2, tmp3, 4*sizeof(su3)/8);
                MPI_Send((void*) tmp2, 4*sizeof(su3)/8, MPI_FLOAT, 0, tag, g_cart_grid);
              }
              else {
                be_to_cpu_assign(tmp, tmp3, 4*sizeof(su3)/8);
                MPI_Send((void*) tmp, 4*sizeof(su3)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
              }
            }
          }
#endif
          tag++;
        }
#ifdef MPI
        MPI_Barrier(g_cart_grid);
#endif
      }
    }
  }

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes, "b");
      fprintf(stdout, "# Time spent writing %s ", measure);
      engineering(measure, tock-tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock-tick), "b/s");
      fprintf(stdout, "# Writing speed: %s", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock-tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
    }
  }
#endif

  return(0);
}
#endif /* HAVE_LIBLEMON */
