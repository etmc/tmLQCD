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
int read_binary_gauge_data(LemonReader * lemonreader, DML_Checksum * checksum, paramsIldgFormat * input, su3 ** const gf)
{
  int t, x, y, z, status = 0;
  int latticeSize[] = {input->lt, input->lx, input->ly, input->lz};
  int scidacMapping[] = {0, 3, 2, 1};
  DML_SiteRank rank;
  MPI_Offset bytes;
  uint64_t fbsu3;
  char * filebuffer = NULL, * current = NULL;
  double tick = 0, tock = 0;
  char measure[64];

  bytes = lemonReaderBytes(lemonreader); /* datalength of ildg-binary-data record in bytes */

  if (bytes != (n_uint64_t)g_nproc * (n_uint64_t)VOLUME * 4 * (n_uint64_t)sizeof(su3) / (input->prec==64 ? 1 : 2)) {
    fprintf(stderr, "Lattice size and precision found in data file do not match those requested at input.\n");
    fprintf(stderr, "Expected LX = %d, LY = %d, LZ = %d, LT = %d, and %s precision.\n", input->lx, input->ly, input->lz, input->lt, (input->prec==64 ? "double" : "single"));
    fprintf(stderr, "Expected %lu bytes, found %lu bytes in gauge file.\n", (unsigned long)(n_uint64_t)g_nproc * (n_uint64_t)VOLUME * 4 * (n_uint64_t)sizeof(su3) / (input->prec==64 ? 1 : 2), (unsigned long)bytes);
    fprintf(stderr, "Check input parameters T, L (LX, LY, LZ) and GaugeConfigReadPrecision.\n");
    return(-3);
  }

  DML_checksum_init(checksum);

  fbsu3 = sizeof(su3);
  if (input->prec == 32) {
    fbsu3 /= 2;
  }
  bytes = 4 * fbsu3;


  if((void*)(filebuffer = malloc(VOLUME * bytes)) == NULL) {
    fprintf (stderr, "malloc errno %d in read_binary_gauge_data, returning without reading gauge file.\n", errno);
    errno = 0;
    return(-1);
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  status = lemonReadLatticeParallelMapped(lemonreader, filebuffer, bytes, latticeSize, scidacMapping);

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();
  }

  if (status != LEMON_SUCCESS) {
    free(filebuffer);
    fprintf(stderr, "Lemon read error occurred with status = %d, while reading in gauge_read_binary.c!\n", status);
    return(-2);
  }

  if (g_debug_level > 0 && g_cart_id == 0) {
    engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes, "b");
    fprintf(stdout, "# Time spent reading %s ", measure);
    engineering(measure, tock - tick, "s");
    fprintf(stdout, "was %s.\n", measure);
    engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock - tick), "b/s");
    fprintf(stdout, "# Reading speed: %s", measure);
    engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock - tick)), "b/s");
    fprintf(stdout, " (%s per MPI process).\n", measure);
    fflush(stdout);
  }


  for (t = 0; t < T; t++) {
    for (z = 0; z < LZ; z++) {
      for (y = 0; y < LY; y++) {
        for (x = 0; x < LX; x++) {
          rank = (DML_SiteRank)(g_proc_coords[1] * LX +
                                (((g_proc_coords[0] * T + t) * g_nproc_z * LZ + g_proc_coords[3] * LZ + z) * g_nproc_y * LY
                                 + g_proc_coords[2] * LY + y) * ((DML_SiteRank)LX * g_nproc_x) + x);
          current = filebuffer + bytes * (x + (y + (t * LZ + z) * LY) * LX);
          DML_checksum_accum(checksum, rank, current, bytes);
          if (input->prec == 32) {
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][1], current            , sizeof(su3) / 8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][2], current +     fbsu3, sizeof(su3) / 8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][3], current + 2 * fbsu3, sizeof(su3) / 8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][0], current + 3 * fbsu3, sizeof(su3) / 8);
          }
          else {
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][1], current            , sizeof(su3) / 8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][2], current +     fbsu3, sizeof(su3) / 8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][3], current + 2 * fbsu3, sizeof(su3) / 8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][0], current + 3 * fbsu3, sizeof(su3) / 8);
          }
        }
      }
    }
  }
  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);
  free(filebuffer);
  return(0);
}
#else /* HAVE_LIBLEMON */
int read_binary_gauge_data(LimeReader * limereader, DML_Checksum * checksum, paramsIldgFormat * input, su3 ** const gf) {

  int t, x, y , z, status=0;
  int latticeSize[] = {input->lt, input->lx, input->ly, input->lz};
  n_uint64_t bytes;
  su3 tmp[4];
  float tmp2[72];
#ifdef MPI
  double tick = 0, tock = 0;
#endif
  char measure[64];
  DML_SiteRank rank;
  DML_checksum_init(checksum);

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
#endif

  bytes = limeReaderBytes(limereader); /* datalength of ildg-binary-data record in bytes */
  if (bytes != (n_uint64_t)g_nproc * (n_uint64_t)VOLUME * 4 * (n_uint64_t)sizeof(su3) / (input->prec==64 ? 1 : 2)) {
    fprintf(stderr, "Lattice size and precision found in data file do not match those requested at input.\n");
    fprintf(stderr, "Expected LX = %d, LY = %d, LZ = %d, LT = %d, and %s precision.\n", input->lx, input->ly, input->lz, input->lt, (input->prec==64 ? "double" : "single"));
    fprintf(stderr, "Expected %lu bytes, found %lu bytes.\n", (unsigned long)(n_uint64_t)g_nproc * (n_uint64_t)VOLUME * 4 * (n_uint64_t)sizeof(su3) / (input->prec==64 ? 1 : 2), (unsigned long)bytes);
    fprintf(stderr, "Check input parameters T, L (LX, LY, LZ) and GaugeConfigReadPrecision.\n");
    return(-3);
  }

  if(input->prec == 32) bytes = (n_uint64_t)2*sizeof(su3);
  else bytes = (n_uint64_t)4*sizeof(su3);
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
#ifdef MPI
        limeReaderSeek(limereader,(n_uint64_t)
                       (((n_uint64_t) g_proc_coords[1]*LX) +
                        ((n_uint64_t) (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
                         + g_proc_coords[2]*LY+y)*LX*g_nproc_x))*bytes,
                       SEEK_SET);
#endif
        for(x = 0; x < LX; x++) {
          rank = (DML_SiteRank) (g_proc_coords[1]*LX +
                                 (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
                                  + g_proc_coords[2]*LY+y)*((DML_SiteRank)LX*g_nproc_x) + x);
          if(input->prec == 32) {
            status = limeReaderReadData(tmp2, &bytes, limereader);
            DML_checksum_accum(checksum, rank, (char *) tmp2, bytes);
          }
          else {
            status = limeReaderReadData(tmp, &bytes, limereader);
            DML_checksum_accum(checksum, rank, (char *) tmp, bytes);
          }
          if(status < 0 && status != LIME_EOR) {
            fprintf(stderr, "LIME read error occurred with status = %d while reading in gauge_read_binary.c!\n", status);
#ifdef MPI
              MPI_Abort(MPI_COMM_WORLD, 1);
              MPI_Finalize();
#endif
            return(-2);
          }
          if(input->prec == 32) {
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][0], &tmp2[3*18], sizeof(su3)/8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][1], &tmp2[0*18], sizeof(su3)/8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][2], &tmp2[1*18], sizeof(su3)/8);
            be_to_cpu_assign_single2double(&gf[ g_ipt[t][x][y][z] ][3], &tmp2[2*18], sizeof(su3)/8);
          }
          else {
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3)/8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3)/8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3)/8);
            be_to_cpu_assign(&gf[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3)/8);
          }
        }
      }
    }
  }

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes, "b");
      fprintf(stdout, "# Time spent reading %s ", measure);
      engineering(measure, tock-tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (tock-tick), "b/s");
      fprintf(stdout, "# Reading speed: %s", measure);
      engineering(measure, latticeSize[0] * latticeSize[1] * latticeSize[2] * latticeSize[3] * bytes / (g_nproc * (tock-tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
    }
  }

  DML_checksum_combine(checksum);
#endif
  return(0);
}
#endif /* HAVE_LIBLEMON */
