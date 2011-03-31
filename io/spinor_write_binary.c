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

#ifdef HAVE_LIBLEMON
int write_binary_spinor_data(spinor * const s, spinor * const r,
                             LemonWriter * lemonwriter, DML_Checksum *checksum, int const prec)
{
  int x, y, z, t, i = 0, xG, yG, zG, tG, status = 0;
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  int scidacMapping[] = {0, 3, 2, 1};
  unsigned long bufoffset = 0;
  char *filebuffer = NULL;
  uint64_t bytes;
  DML_SiteRank rank;
  double tick = 0, tock = 0;
  char measure[64];
  spinor *p = NULL;

  DML_checksum_init(checksum);
  bytes = (uint64_t)sizeof(spinor);
  if (prec == 32) {
    bytes /= 2;
  }
  if((void*)(filebuffer = malloc(VOLUME * bytes)) == NULL) {
    fprintf (stderr, "malloc errno in write_binary_spinor_data_parallel: %d\n", errno);
    fflush(stderr);
    errno = 0;
    /* do we need to abort here? */
    return 1;
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
          i = g_lexic2eosub[g_ipt[t][x][y][z]];
          if ((z  + zG + y  + yG +
               x  + xG + t + tG) % 2 == 0)
            p = s;
          else
            p = r;

          if (prec == 32)
            be_to_cpu_assign_double2single((float*)(filebuffer + bufoffset), (double*)(p + i), sizeof(spinor) / 8);
          else
            be_to_cpu_assign((double*)(filebuffer + bufoffset), (double*)(p + i),  sizeof(spinor) / 8);
            DML_checksum_accum(checksum, rank, (char*) filebuffer + bufoffset, bytes);
            bufoffset += bytes;
        }
      }
    }
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  status = lemonWriteLatticeParallelMapped(lemonwriter, filebuffer, bytes, latticeSize, scidacMapping);

  if (status != LEMON_SUCCESS)
  {
    free(filebuffer);
    fprintf(stderr, "LEMON write error occurred with status = %d, while in write_binary_spinor_data (spinor_write_binary.c)!\n", status);
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

  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);

  free(filebuffer);
  return 0;

}

#else /* HAVE_LIBLEMON */
int write_binary_spinor_data(spinor * const s, spinor * const r, LimeWriter * limewriter, DML_Checksum * checksum, const int prec)
{
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0, status=0;
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  int coords[4];
  n_uint64_t bytes;
  DML_SiteRank rank;
#ifdef MPI
  double tick = 0, tock = 0;
  char measure[64];
  MPI_Status mstatus;
#endif
  DML_checksum_init(checksum);

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
#endif

  if(prec == 32) bytes = (n_uint64_t)sizeof(spinor)/2;
  else bytes = (n_uint64_t)sizeof(spinor);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    t = t0 - T*g_proc_coords[0];
    coords[0] = t0 / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
        Y = y - g_proc_coords[2]*LY;
        coords[2] = y / LY;
        for(x = 0; x < LX*g_nproc_x; x++) {
          X = x - g_proc_coords[1]*LX;
          coords[1] = x / LX;
#ifdef MPI
          MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
          if(g_cart_id == id) {
            i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
            if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY
                + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
              p = s;
            }
            else {
              p = r;
            }
          }
          if(g_cart_id == 0) {
            /* Rank should be computed by proc 0 only */
            rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);

            if(g_cart_id == id) {
              if(prec == 32) {
                be_to_cpu_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8);
                DML_checksum_accum(checksum,rank,(char *) tmp2,sizeof(spinor)/2);
                status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
              }
              else {
                be_to_cpu_assign(tmp, p + i , sizeof(spinor)/8);
                DML_checksum_accum(checksum,rank,(char *) tmp,sizeof(spinor));
                status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
              }
            }
#ifdef MPI
            else{
              if(prec == 32) {
                MPI_Recv((void*)tmp2, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &mstatus);
                DML_checksum_accum(checksum,rank,(char *) tmp2, sizeof(spinor)/2);
                status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
              }
              else {
                MPI_Recv((void*)tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mstatus);
                DML_checksum_accum(checksum,rank,(char *) tmp, sizeof(spinor));
                status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
              }
            }
#endif
            if(status < 0 ) {
              fprintf(stderr, "LIME write error occurred with status = %d, while in write_binary_spinor_data (spinor_write_binary.c)!\n", status);
#ifdef MPI
              MPI_Abort(MPI_COMM_WORLD, 1);
              MPI_Finalize();
#endif
              exit(500);
            }
          }
#ifdef MPI
          else{
            if(g_cart_id == id){
              if(prec == 32) {
                be_to_cpu_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8);
                MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
              }
              else {
                be_to_cpu_assign(tmp, p + i, sizeof(spinor)/8);
                MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
              }
            }
          }
#endif
          tag++;
        }
#ifdef MPI
        MPI_Barrier(g_cart_grid);
#endif
        tag=0;
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
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock - tick), "b/s");
      fprintf(stdout, "# Writing speed: %s", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }
#endif
  return(0);
}
#endif /* HAVE_LIBLEMON */



#ifdef HAVE_LIBLEMON
int write_binary_spinor_data_l(spinor * const s,
                             LemonWriter * lemonwriter, DML_Checksum *checksum, int const prec)
{
  int x, y, z, t, i = 0, xG, yG, zG, tG, status = 0;
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  int scidacMapping[] = {0, 3, 2, 1};
  unsigned long bufoffset = 0;
  char *filebuffer = NULL;
  uint64_t bytes;
  DML_SiteRank rank;
  double tick = 0, tock = 0;
  char measure[64];

  DML_checksum_init(checksum);
  bytes = (uint64_t)sizeof(spinor);
  if (prec == 32) {
    bytes /= 2;
  }
  if((void*)(filebuffer = malloc(VOLUME * bytes)) == NULL) {
    fprintf (stderr, "malloc errno in write_binary_spinor_data_parallel: %d\n", errno);
    fflush(stderr);
    errno = 0;
    /* do we need to abort here? */
    return 1;
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
          i = g_ipt[t][x][y][z];

          if (prec == 32)
            be_to_cpu_assign_double2single((float*)(filebuffer + bufoffset), (double*)(s + i), sizeof(spinor) / 8);
          else
            be_to_cpu_assign((double*)(filebuffer + bufoffset), (double*)(s + i),  sizeof(spinor) / 8);
            DML_checksum_accum(checksum, rank, (char*) filebuffer + bufoffset, bytes);
            bufoffset += bytes;
        }
      }
    }
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  status = lemonWriteLatticeParallelMapped(lemonwriter, filebuffer, bytes, latticeSize, scidacMapping);

  if (status != LEMON_SUCCESS)
  {
    free(filebuffer);
    fprintf(stderr, "LEMON write error occurred with status = %d, while in write_binary_spinor_data_l (spinor_write_binary.c)!\n", status);
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

  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);

  free(filebuffer);
  return 0;

}

#else /* HAVE_LIBLEMON */
int write_binary_spinor_data_l(spinor * const s, LimeWriter * limewriter, DML_Checksum * checksum, const int prec)
{
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0, status=0;
  int latticeSize[] = {T_global, g_nproc_x*LX, g_nproc_y*LY, g_nproc_z*LZ};
  spinor tmp[1];
  float tmp2[24];
  int coords[4];
  n_uint64_t bytes;
  DML_SiteRank rank;
#ifdef MPI
  double tick = 0, tock = 0;
  char measure[64];
  MPI_Status mstatus;
#endif
  DML_checksum_init(checksum);

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
#endif

  if(prec == 32) bytes = (n_uint64_t)sizeof(spinor)/2;
  else bytes = (n_uint64_t)sizeof(spinor);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    t = t0 - T*g_proc_coords[0];
    coords[0] = t0 / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
        Y = y - g_proc_coords[2]*LY;
        coords[2] = y / LY;
        for(x = 0; x < LX*g_nproc_x; x++) {
          X = x - g_proc_coords[1]*LX;
          coords[1] = x / LX;
#ifdef MPI
          MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
          if(g_cart_id == id) {
            i = g_ipt[t][X][Y][Z];
          }
          if(g_cart_id == 0) {
            /* Rank should be computed by proc 0 only */
            rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);

            if(g_cart_id == id) {
              if(prec == 32) {
                be_to_cpu_assign_double2single((float*)tmp2, s + i, sizeof(spinor)/8);
                DML_checksum_accum(checksum,rank,(char *) tmp2,sizeof(spinor)/2);
                status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
              }
              else {
                be_to_cpu_assign(tmp, s + i , sizeof(spinor)/8);
                DML_checksum_accum(checksum,rank,(char *) tmp,sizeof(spinor));
                status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
              }
            }
#ifdef MPI
            else{
              if(prec == 32) {
                MPI_Recv((void*)tmp2, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &mstatus);
                DML_checksum_accum(checksum,rank,(char *) tmp2, sizeof(spinor)/2);
                status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
              }
              else {
                MPI_Recv((void*)tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mstatus);
                DML_checksum_accum(checksum,rank,(char *) tmp, sizeof(spinor));
                status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
              }
            }
#endif
            if(status < 0 ) {
              fprintf(stderr, "LIME write error occurred with status = %d, while in write_binary_spinor_data_l (spinor_write_binary.c)!\n", status);
#ifdef MPI
              MPI_Abort(MPI_COMM_WORLD, 1);
              MPI_Finalize();
#endif
              exit(500);
            }

          }
#ifdef MPI
          else{
            if(g_cart_id == id){
              if(prec == 32) {
                be_to_cpu_assign_double2single((float*)tmp2, s + i, sizeof(spinor)/8);
                MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
              }
              else {
                be_to_cpu_assign(tmp, s + i, sizeof(spinor)/8);
                MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
              }
            }
          }
#endif
          tag++;
        }
#ifdef MPI
        MPI_Barrier(g_cart_grid);
#endif
        tag=0;
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
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock - tick), "b/s");
      fprintf(stdout, "# Writing speed: %s", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }
#endif
  return(0);
}
#endif /* HAVE_LIBLEMON */
