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
/* $Id$ */

/****************************************************
 * IO routines:
 *
 * write_lime_gauge_field_old(char * filename, const double plaq, const int counter)
 *   write gauge field in ILDG format
 *
 * read_lime_gauge_field_old(char * filename)
 *   read gauge field in ILDG format
 *
 * write_spinorfield_eo_time_p(spinor * const s, spinor * const r, char * filename, const int append)
 *   write a full spinor from one even and one odd spinorfield
 *
 * read_spinorfield_eo_time(spinor * const s, spinor * const r, char * filename)
 *   read a full spinor into one even and one odd spinorfields
 *
 * write_gauge_field(char * filename)
 *   writes gauge field configuration to file
 *   with name filename.
 *
 * read_gauge_field(char * filename)
 *   reads gauge field configuration from file
 *   with name filename.
 *
 * int read_eospinor(spinor * const s, char * filename)
 *   Read an even or odd spinor from file (eigenvector)
 *
 * int write_eospinor(spinor * const s, char * filename, const double evalue, const double prec, const int nstore)
 *   Read an even or odd spinor from file (eigenvector)
 *
 * Autor: Ines Wetzorke <wetzorke@ifh.de>
 *        Carsten Urbach <urbach@ifh.de>
 *
 * Addition of isnan in function read_spinorfield_cm_single to detect a
 * wrong read order and creation of read_spinorfield_cm_swap_single to
 * swap while reading if necessary.
 *
 * Author: Remi Baron <remi.baron@cea.fr> Jul 2007
 *
 ****************************************************/

/*
 * Note:
 * Required version of lime: >= 1.2.3
 * n_uint64_t is a lime defined type!!
 *
 */

#define _FILE_OFFSET_BITS 64

#include"lime.h" 
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<sys/time.h> 
#include<sys/types.h>
#ifdef MPI
# include<unistd.h> 
#endif
#include<math.h>
#include"global.h"
#include"su3.h"
#include"lime.h" 
#include"read_input.h"
#include<io/utils.h>
#include"io.h"

#define MAXBUF 1048576



/*************************************************
 *
 * This routine writes an even or odd spinor-field
 * so really of size VOLUME/2
 *
 * used for instance for storing eigenvectors
 * of the precoditioned matrix
 *
 *************************************************/

int write_eospinor(spinor * const s, char * filename, 
		   const double evalue, const double prec, const int nstore) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0;
  int ME_flag=0, MB_flag=0, status=0;
  spinor tmp[1];
  int coords[4];
  char message[500];
  n_uint64_t bytes;
#ifdef MPI
  MPI_Status mpistatus;
#endif

  if(g_cart_id == 0){  
    if(g_kappa > 0. || g_kappa < 0.) {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1, PACKAGE_VERSION);
    }
    else {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu, g_rgi_C1, PACKAGE_VERSION);
    }
    bytes = strlen( message );
    
    if((ofs = fopen(filename, "w")) == (FILE*)NULL) {
      fprintf(stderr, "Error writing eigenvector to file %s!\n", filename);
      return(-1);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }

    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (xlf-info) error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);

    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2;
    MB_flag=0; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "eospinor-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (eospinor-binary-data) error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX*g_nproc_x; x++){
    X = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++){
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++){
	Z = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++){
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	  if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    if(g_cart_id == 0) {
	      if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
		byte_swap_assign(tmp, s + i , sizeof(spinor)/8);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
#else
		status = limeWriteRecordData((void*)(s+i), &bytes, limewriter);
#endif
	      }
#ifdef MPI
	      else {
		MPI_Recv(tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpistatus);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#endif
	      if(status < 0 ) {
		fprintf(stderr, "LIME write error %d\n", status);
#ifdef MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
#endif
		exit(500);
	      }
	    }
#ifdef MPI
	    else {
	      if(g_cart_id == id) {
#  ifndef WORDS_BIGENDIAN
		byte_swap_assign(tmp, s + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  else
		MPI_Send((void*) (s + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
#  endif		  
	      }
	    }
#endif
	    tag++;
	  }
	}
#ifdef MPI
 	MPI_Barrier(g_cart_grid); 
#endif
	tag=0;
      }
    }
  }
  if(g_cart_id == 0) {
    if(ferror(ofs)) {
      fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
    }
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  return(0);
}

int read_eospinor(spinor * const s, char * filename) {
  FILE * ifs;
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
#ifdef MPI
  int position;
#endif
#ifndef WORDS_BIGENDIAN
  spinor tmp[1];
#endif
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Error opening file %s\n", filename);
    }
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Unable to open LimeReader\n");
    }
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("eospinor-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no eospinor-binary-data record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2) {
    if(g_proc_id == 0) {
      fprintf(stderr, "wrong length in eospinor: %d. Aborting read!\n", (int)bytes);
    }
    return(-1);
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	limeReaderSeek(limereader, (n_uint64_t)
		       (g_proc_coords[0]*T+
			(((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
			 + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
		       SEEK_SET);
#endif
	for(t = 0; t < T; t++){
	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[3]*LZ+g_proc_coords[2]*LY
	      +g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    
#ifndef WORDS_BIGENDIAN
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(s + i, tmp, sizeof(spinor)/8);
#else
	    status = limeReaderReadData((s+i), &bytes, limereader);
#endif
	    if(status < 0 && status != LIME_EOR) {
	      fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", status, filename);
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
#endif
	      exit(500);
	    }
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}




int read_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename, 
			       const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    return(-1);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	fseek(ifs,
	      (g_proc_coords[0]*T+
	       (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
		+ g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
	      SEEK_SET);
#endif
	for(t = 0; t < T; t++) {

	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[0]*T+g_proc_coords[1]*LX+
	      g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  
	  if(ts == t || ts < 0 || ts >= T){
	    /* Read the data */
	    fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Test if we read the data with the correct endian order */
            if(isnan(tmp[0]) || isnan(tmp[1]) || isnan(tmp[2]) || isnan(tmp[3]) || isnan(tmp[4]) || isnan(tmp[5]) ||
            isnan(tmp[6]) || isnan(tmp[7]) || isnan(tmp[8]) || isnan(tmp[9]) || isnan(tmp[10]) || isnan(tmp[11]) ||
            isnan(tmp[12]) || isnan(tmp[13]) || isnan(tmp[14]) || isnan(tmp[15]) || isnan(tmp[16]) || isnan(tmp[17]) ||
            isnan(tmp[18]) || isnan(tmp[19]) || isnan(tmp[20]) || isnan(tmp[21]) || isnan(tmp[22]) || isnan(tmp[23]))
            {
              if(g_proc_id == 0)
              {
                if(big_endian())
                  printf("\nBig endian order gives some NaN. Trying little endian order instead...\n\n");
                else
                  printf("\nLittle endian order gives some NaN. Trying big endian order instead...\n\n");
              }

              fclose(ifs);
              return read_spinorfield_cm_swap_single(s,r,filename,ts,vol);
            }
	    single2double_cm(p+i, tmp);
	  }
	  else {
	    if(vol > 0) {
	      fread(tmp, sizeof(spinor)/2, 1, ifs);
	    }
	    /* Padding with zeros */
	    zero_spinor(p+i);
	  }
	}
      }
    }
  }
  fclose(ifs);
  return(0);
}

int read_spinorfield_cm_swap_single(spinor * const s, spinor * const r, char * filename,
                               const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
        fseek(ifs,
              (g_proc_coords[0]*T+
               (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
                + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
              SEEK_SET);
#endif
        for(t = 0; t < T; t++) {

          i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
          if((t+x+y+z+
              g_proc_coords[0]*T+g_proc_coords[1]*LX+
              g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
            p = s;
          }
          else {
            p = r;
          }

          if(ts == t || ts < 0 || ts >= T){
            /* Read the data */
            fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Swap and convert from single to double precision */
            byte_swap_assign_single2double(p+i, tmp, sizeof(spinor)/8);
          }
          else {
            if(vol > 0) {
              fread(tmp, sizeof(spinor)/2, 1, ifs);
            }
            /* Padding with zeros */
            zero_spinor(p+i);
          }
        }
      }
    }
  }
  fclose(ifs);
  return(0);
}


int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename) {

  FILE * ofs = NULL;
  int t, x, y , z, i = 0;
  int t0, X, Y, Z, id = 0;
  spinor * p = NULL;
  float tmp[24];
  int coords[4];
#ifdef MPI
  int  tag = 0;
  MPI_Status status;
#endif
  
  if(g_cart_id == 0) {
    ofs = fopen(filename, "w");
    printf("# Writing in cmi format (32 Bit) to file %s\n", filename);
  }

  for(x = 0; x < LX*g_nproc_x; x++) {
    X = x - LX*g_proc_coords[1];
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      Y = y - LY*g_proc_coords[2];
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	Z = z - LZ*g_proc_coords[3];
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++) {
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
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
	  if(g_cart_id == 0){
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	    }
#ifdef MPI
	    else {
	      MPI_Recv(tmp, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &status);
	    }
#endif
#ifndef WORDS_BIGENDIAN
            byte_swap(tmp,sizeof(spinor)/8);
#endif
	    fwrite(tmp, sizeof(float), 24, ofs);
	  }
#ifdef MPI
	  else {
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	      MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	    }
	  }
	  tag++;
#endif
	}
#ifdef MPI
	MPI_Barrier(g_cart_grid); 
	tag=0;
#endif
      }
    }
  }
  if(g_cart_id == 0) {
    fclose(ofs);
  }
  return(0);
}




