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

#define _FILE_OFFSET_BITS 64

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_LIBLEMON
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>
#include "global.h"
#include "su3.h"

#include "parallel_io.h"

#include "read_input.h"
#include "io_utils.h"

#define PACKAGE_VERSION "5.0.2"

#include "ranlxd.h"

int read_lemon_gauge_field_parallel(char *filename, int *ranlxd_state)
{
  MPI_File ifs;
  int status;
  char * header_type;
  LemonReader * lemonreader;
  DML_Checksum checksum_read;
  DML_Checksum checksum_calc;
  int DML_read_flag = 0;
  int gauge_read_flag = 0;
  int ranlxd_read_flag = 0;

  status = MPI_File_open(g_cart_grid, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &ifs);
  if (status)
  {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonreader = lemonCreateReader(&ifs, g_cart_grid);
  if( lemonreader == (LemonReader *)NULL )
  {
    fprintf(stderr, "Unable to open LemonReader\n");
    MPI_File_close(&ifs);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  /* This is added for convenience - not sure about calling conventions here. */
  while( (status = lemonReaderNextRecord(lemonreader)) == LEMON_SUCCESS )
  {
    header_type = lemonReaderType(lemonreader);
    if (strcmp("ranluxd-state", header_type) == 0)
    {
      read_rlxd_state_parallel(lemonreader, ranlxd_state);
      ranlxd_read_flag = 1;
    }
    else if (strcmp("ildg-binary-data", header_type) == 0)
    {
      read_binary_gauge_data_parallel(lemonreader, &checksum_calc);
      gauge_read_flag = 1;
    }
    else if (strcmp("scidac-checksum", header_type) == 0)
    {
      if (g_debug_level > 1)
      {
        read_checksum_parallel(lemonreader, &checksum_read);
        DML_read_flag = 1;
      }
    }
    lemonReaderCloseRecord(lemonreader);
  }

  if (!gauge_read_flag)
  {
    if (g_cart_id == 0)
    {
      fprintf(stderr, "Did not find gauge record in %s.\n", filename);
      fprintf(stderr, "Panic! Aborting...\n");
      fflush(stderr);
    }
    MPI_File_close(lemonreader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(501);
  }

  if (!ranlxd_read_flag)
  {
    if (g_cart_id == 0)
    {
      fprintf(stderr, "Did not find a ranluxd-state record in %s.\n", filename);
      fprintf(stderr, "Leaving state uninitialized.\n");
      fflush(stderr);
    }
  }

  if (g_debug_level > 1 && g_cart_id == 0)
  {
    printf("# checksum for gaugefield %s\n", filename);
    printf("# calculated: %#x %#x.\n", checksum_calc.suma, checksum_calc.sumb);
    if (DML_read_flag)
      printf("# read:       %#x %#x.\n", checksum_read.suma, checksum_read.sumb);
    else
      printf("# A DML checksum record was not found in %s.\n", filename);
    fflush(stderr);
  }

  MPI_File_close(&ifs);
  return 0; /*TODO Error handling, etc.*/
}

int read_binary_gauge_data_parallel(LemonReader * lemonreader, DML_Checksum * checksum)
{
  int t, x, y, z, status=0;
  int latticeSize[] = {L, L, L, T_global};
  int prec;
  DML_SiteRank rank;
  MPI_Offset bytes;
  uint64_t fbsu3;
  char *filebuffer, *current;

  bytes = lemonReaderBytes(lemonreader);

  if (bytes == g_nproc * VOLUME * 4 * sizeof(su3))
    prec = 64;
  else if (bytes == g_nproc * VOLUME * 4 * sizeof(su3) / 2)
    prec = 32;
  else
  {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu).s\n", (unsigned long)bytes);
    fprintf(stderr, "Panic! Aborting...\n");
    fflush( stdout );
    MPI_File_close(lemonreader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(501);
  }

  if(g_cart_id == 0 && g_debug_level > 2)
    printf("# %d Bit precision read.\n", prec);

  DML_checksum_init(checksum);

  fbsu3 = sizeof(su3);
  if (prec == 32)
    fbsu3 /= 2;
  bytes = 4 * fbsu3;

  filebuffer = malloc(VOLUME * bytes);

  lemonReadLatticeParallel(lemonreader, filebuffer, bytes, latticeSize);

  if(status < 0 && status != LEMON_EOR)
  {
    fprintf(stderr, "LEMON read error occured with status = %d while reading!\nPanic! Aborting...\n", status);
    MPI_File_close(lemonreader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  for(t = 0; t < T; t++)
    for(z = 0; z < LZ; z++)
      for(y = 0; y < LY; y++)
        for(x = 0; x < LX; x++)
        {
          rank = (DML_SiteRank) (g_proc_coords[1]*LX +
              (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
                + g_proc_coords[2]*LY+y)*((DML_SiteRank)LX*g_nproc_x) + x);
          current = filebuffer + bytes * (x + (y + (t * LZ + z) * LY) * LX);
          DML_checksum_accum(checksum, rank, current, bytes);
#ifndef WORDS_BIGENDIAN
          if(prec == 32)
          {
            byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], current            , sizeof(su3)/8);
            byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], current +     fbsu3, sizeof(su3)/8);
            byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], current + 2 * fbsu3, sizeof(su3)/8);
            byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], current + 3 * fbsu3, sizeof(su3)/8);
          }
          else
          {
            /* NOTE Changed size of element here - seemed wrong. Check!!! */
            byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][1], current            , sizeof(su3)/4);
            byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][2], current +     fbsu3, sizeof(su3)/4);
            byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][3], current + 2 * fbsu3, sizeof(su3)/4);
            byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][0], current + 3 * fbsu3, sizeof(su3)/4);
          }
#else
          if(prec == 32)
          {
            single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], current            , sizeof(su3)/8);
            single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], current +     fbsu3, sizeof(su3)/8);
            single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], current + 2 * fbsu3, sizeof(su3)/8);
            single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], current + 3 * fbsu3, sizeof(su3)/8);
          }
          else
          {
            memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][0], current            , sizeof(su3));
            memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][1], current +     fbsu3, sizeof(su3));
            memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][2], current + 2 * fbsu3, sizeof(su3));
            memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][3], current + 3 * fbsu3, sizeof(su3));
          }
#endif
        }
  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);
  return(0);
}

int read_checksum_parallel(LemonReader * lemonreader, DML_Checksum * checksum)
{
  char *message;
  uint64_t bytes, bytes_read;

  bytes = lemonReaderBytes(lemonreader);
  bytes_read = bytes;
  message = (char*)malloc(bytes + 1); /* Allow for an explicit closing byte. */

  lemonReaderReadData(message, &bytes_read, lemonreader);
  message[bytes] = '\0';

  /* So, yeah, that means we're going to have to parse the XML */
  sscanf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                   "<scidacChecksum>\n"
                   "  <version>1.0</version>\n"
                   "  <suma>%#x</suma>\n"
                   "  <sumb>%#x</sumb>\n"
                   "</scidacChecksum>", &checksum->suma, &checksum->sumb);

  return(0);
}

int read_rlxd_state_parallel(LemonReader *lemonreader, int *state)
{
  /* We assume we are passed a reader set to the appropriate record */
  int status;
  uint64_t len;

  len = lemonReaderBytes(lemonreader);
  if((int)len != rlxdsize*sizeof(int))
  {
    fprintf(stderr, "Wrong size for ranluxd-state (len=%d!=%d)\n", (int)len, (int)(rlxdsize*sizeof(int)));
    fprintf(stderr, "Panic! Aborting...!\n");
    fflush( stdout );
    MPI_File_close(lemonreader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(501);
  }
  status = lemonReaderReadData(state, &len, lemonreader);

  if(status < 0 && status != LEMON_EOR)
  {
    fprintf(stderr, "LEMON read error occured with status = %d!\nPanic! Aborting...\n", status);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  #ifndef WORDS_BIGENDIAN
  byte_swap(state, rlxdsize);
  #endif

  return(0);
}

int write_lemon_gauge_field_parallel(char * filename, const double plaq, const int counter, const int prec) {
  MPI_File ofs;
  LemonWriter * lemonwriter = NULL;
  LemonRecordHeader * lemonheader = NULL;

  int status=0;
  uint64_t bytes;
  DML_Checksum checksum;

  MPI_File_open(g_cart_grid, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ofs);

  lemonwriter = lemonCreateWriter( &ofs, g_cart_grid );
  if(lemonwriter == (LemonWriter*)NULL) {
    fprintf(stderr, "LEMON error in file %s for writing!\nPanic! Aborting...\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  write_xlf_info_parallel(lemonwriter, plaq, counter);

  bytes = ((uint64_t)LX*g_nproc_x)*((uint64_t)LY*g_nproc_y)*((uint64_t)LZ*g_nproc_z)*((uint64_t)T*g_nproc_t)*((uint64_t)4*sizeof(su3));
  if(prec == 32) bytes = bytes/((uint64_t)2);

  write_ildg_format_parallel(lemonwriter, prec);

  lemonheader = lemonCreateHeader(0, 0, "ildg-binary-data", bytes);
  status = lemonWriteRecordHeader(lemonheader, lemonwriter);
  if(status < 0 ) {
    fprintf(stderr, "LEMON write header error %d\n", status); fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader(lemonheader);

  write_binary_gauge_data_parallel(lemonwriter, prec, &checksum);

  if(g_cart_id == 0) {
    fprintf(stderr, "Checksum A: %#x \nChecksum B: %#x\n", checksum.suma, checksum.sumb);
    fflush(stderr);
  }
  write_checksum_parallel(lemonwriter, &checksum);
  write_rlxd_state_parallel(lemonwriter);
  MPI_Barrier(MPI_COMM_WORLD);
  lemonDestroyWriter( lemonwriter );
  MPI_File_close(&ofs);
  return(0);
}

int write_binary_gauge_data_parallel(LemonWriter * lemonwriter, const int prec, DML_Checksum * ans) {
  int x, X, y, Y, z, Z, tt, t0, tag=0, id=0;
  su3 tmp3[4];
  int coords[4];

  int globaldims[] = {L, L, L, T_global};
  unsigned long bufoffset;
  char *filebuffer;
  uint64_t bytes;
  DML_SiteRank rank;
  DML_checksum_init(ans);
  bytes = (uint64_t)sizeof(su3) * (prec == 32 ? 2 :4);
  bufoffset = 0;
  filebuffer = (char*)malloc(bytes * VOLUME);
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
          MPI_Cart_rank(g_cart_grid, coords, &id);
          rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);
          if(g_cart_id == id) {
            memcpy(&tmp3[0], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][1], sizeof(su3));
            memcpy(&tmp3[1], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][2], sizeof(su3));
            memcpy(&tmp3[2], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][3], sizeof(su3));
            memcpy(&tmp3[3], &g_gauge_field[ g_ipt[tt][X][Y][Z] ][0], sizeof(su3));
            #ifndef WORDS_BIGENDIAN
            if(prec == 32)
              byte_swap_assign_double2single(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            else
              byte_swap_assign(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            #else
            if(prec == 32)
              double2single(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            else
              memcpy(filebuffer + bufoffset, tmp3, 4*sizeof(su3)/8);
            #endif
            DML_checksum_accum(ans, rank, (char*) filebuffer + bufoffset, bytes);

            bufoffset += bytes;
          }
          tag++;
        }
      }
    }
  }
  lemonWriteLatticeParallel(lemonwriter, filebuffer, bytes, globaldims);
  lemonWriterCloseRecord(lemonwriter);

  DML_global_xor(&ans->suma);
  DML_global_xor(&ans->sumb);
  return(0);
}

int write_ildg_format_parallel(LemonWriter *writer, const int prec)
{
  uint64_t bytes;
  int status = 0;
  LemonRecordHeader *header;
  char *buf;

  buf = (char*)malloc(512);
  sprintf(buf, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n"
               "            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
               "            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n"
               "  <version> 1.0 </version>\n"
               "  <field> su3gauge </field>\n"
               "  <precision> %d </precision>\n"
               "  <lx> %d </lx>\n"
               "  <ly> %d </ly>\n"
               "  <lz> %d </lz>\n"
               "  <lt> %d </lt>\n"
               "</ildgFormat>", prec, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z, T*g_nproc_t);

  bytes = strlen(buf);
  header = lemonCreateHeader(1, 0, "ildg-format", bytes);
  if(header == (LemonRecordHeader*)NULL)
  {
    fprintf(stderr, "LEMON create header ildg-format error\nPanic! Aborting...\n");
    exit(500);
  }
  status = lemonWriteRecordHeader(header, writer);
  if(status < 0 )
  {
    fprintf(stderr, "LEMON write header ildg-format error %d\nPanic! Aborting...\n", status);
    exit(500);
  }
  lemonDestroyHeader(header);
  lemonWriteRecordData(buf, &bytes, writer);
  lemonWriterCloseRecord(writer);
  free(buf);
  return 0;
}

int write_checksum_parallel(LemonWriter * lemonwriter, DML_Checksum * checksum)
{
  LemonRecordHeader * lemonheader = NULL;
  int status = 0;
  int MB_flag = 0, ME_flag = 1;
  char message[500];
  uint64_t bytes;

  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                   "<scidacChecksum>\n"
                   "  <version>1.0</version>\n"
                   "  <suma>%#x</suma>\n"
                   "  <sumb>%#x</sumb>\n"
                   "</scidacChecksum>", checksum->suma, checksum->sumb);
  bytes = strlen( message );
  lemonheader = lemonCreateHeader(MB_flag, ME_flag, "scidac-checksum", bytes);
  status = lemonWriteRecordHeader( lemonheader, lemonwriter);
  if(status < 0 ) {
    fprintf(stderr, "LEMON write header error %d\n", status);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader( lemonheader );
  lemonWriteRecordData(message, &bytes, lemonwriter);
  lemonWriterCloseRecord(lemonwriter);
  return(0);
}

int write_xlf_info_parallel(LemonWriter * lemonwriter, const double plaq, const int counter){
  LemonRecordHeader * lemonheader = NULL;
  int status=0;
  char message[5000];
  uint64_t bytes;
  struct timeval t1;

  gettimeofday(&t1,NULL);
  if(g_kappa > 0. || g_kappa < 0.) {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n time = %ld\n hmcversion = %s\n mubar = %f\n epsilonbar = %f\n date = %s",
        plaq, counter, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1,t1.tv_sec, PACKAGE_VERSION,
        g_mubar/2./g_kappa, g_epsbar/2./g_kappa, ctime(&t1.tv_sec));
  }
  else {
    sprintf(message,"\n plaquette = %e\n trajectory nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f\n date = %s",
        plaq, counter, g_beta, g_kappa, g_mu, g_rgi_C1, ctime(&t1.tv_sec));
  }
  bytes = strlen( message );
  lemonheader = lemonCreateHeader(1, 1, "xlf-info", bytes);  /* Message end and Message begin flags are 1 */
  status = lemonWriteRecordHeader( lemonheader, lemonwriter);

  if(status < 0 ) {
    fprintf(stderr, "LEMON write header error %d\n", status);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader( lemonheader );
  lemonWriteRecordData(message, &bytes, lemonwriter);
  lemonWriterCloseRecord(lemonwriter);

  return(0);
}

int write_rlxd_state_parallel(LemonWriter *writer)
{
  uint64_t len;
  LemonRecordHeader *header = NULL;
  int * state;
  int status;

  len = rlxdsize * sizeof(int);

  state = (int*)malloc(len);
  if (g_proc_id == 0)
    rlxd_get(state);
  MPI_Barrier(writer->cartesian);
  MPI_Bcast(state, rlxdsize, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Barrier(writer->cartesian);

  header = lemonCreateHeader(1, 1, "ranluxd-state", len);
  status = lemonWriteRecordHeader(header, writer);
  if(status < 0 )
  {
    fprintf(stderr, "LEMON write header error %d\n", status);
    MPI_File_close(writer->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  lemonDestroyHeader(header);

#ifndef WORDS_BIGENDIAN
  byte_swap(state, rlxdsize);
#endif
  status = lemonWriteRecordData(state, &len, writer);
  if(status < 0 )
  {
    fprintf(stderr, "LEMON write error %d\n", status);
    MPI_File_close(writer->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  free(state);
  return(0);
}

#endif
