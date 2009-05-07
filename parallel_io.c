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
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<sys/time.h> 
#include<sys/types.h>
#include<mpi.h>
#include<unistd.h> 
#include<math.h>
#include"global.h"
#include"su3.h"

#include"read_input.h"
#include"io_utils.h"

#include"io.h"
#include"parallel_io.h"


int write_lime_gauge_field_parallel(char * filename, const double plaq, const int counter, const int prec) {
  MPI_File ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;

  int ME_flag=0, MB_flag=0, status=0;
  n_uint64_t bytes;
  DML_Checksum checksum;
  MPI_Info info;

  MPI_Info_create(&info);
  MPI_File_open(g_cart_grid, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &ofs);

  limewriter = limeCreateWriter( &ofs, g_cart_grid );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\nAborting...\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  write_parallel_xlf_info(plaq, counter, limewriter);
  
  bytes = ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3));
  if(prec == 32) bytes = bytes/((n_uint64_t)2);
  MB_flag=1; ME_flag=0;
  limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", bytes);

  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status); fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  limeDestroyHeader( limeheader );

  write_binary_gauge_data_parallel(limewriter, prec, &checksum);

  if(g_cart_id == 0) {
    fprintf(stderr, "Checksum A: %#x \nChecksum B: %#x\n", checksum.suma, checksum.sumb);
    fflush(stderr);
  }
  write_parallel_checksum(limewriter, &checksum);
  MPI_Barrier(MPI_COMM_WORLD);
  limeDestroyWriter( limewriter );
  MPI_File_close(&ofs);
  return(0);
}

int write_binary_gauge_data_parallel(LimeWriter * limewriter, const int prec, DML_Checksum * ans) {
  int x, X, y, Y, z, Z, tt, t0, tag=0, id=0;
  su3 tmp3[4];
  int coords[4];

  int globaldims[] = {L, L, L, T_global};
  unsigned long bufoffset;
  char *filebuffer;
  n_uint64_t bytes;
  DML_SiteRank rank;
  DML_checksum_init(ans);
  bytes = (n_uint64_t)sizeof(su3) * (prec == 32 ? 2 :4);
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
  limeWriteLatticeParallel(limewriter, filebuffer, bytes, globaldims);
  limeWriterCloseRecord(limewriter);

  DML_global_xor(&ans->suma);
  DML_global_xor(&ans->sumb);
  return(0);
}

int write_checksum_parallel(LimeWriter * limewriter, DML_Checksum * checksum){
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int MB_flag = 0, ME_flag = 1;
  char message[500];
  n_uint64_t bytes;

  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum><version>1.0</version><suma>%#x</suma><sumb>%#x</sumb></scidacChecksum>", (*checksum).suma, (*checksum).sumb);
  bytes = strlen( message );
  limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-checksum", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);
  limeWriterCloseRecord(limewriter);
  return(0);
}

int write_xlf_info_parallel(const double plaq, const int counter, LimeWriter * limewriter){
  LimeRecordHeader * limeheader = NULL;
  int status=0;
  char message[5000];
  n_uint64_t bytes;
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
  limeheader = limeCreateHeader(1, 1, "xlf-info", bytes);  /* Message end and Message begin flags are 1 */
  status = limeWriteRecordHeader( limeheader, limewriter);

  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);
  limeWriterCloseRecord(limewriter);

  return(0);
}

