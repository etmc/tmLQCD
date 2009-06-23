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
#include"dml.h"
#include"io_utils.h"

int isnan_f  (float       x) { return x != x; }
int isnan_d  (double      x) { return x != x; }
int isnan_ld (long double x) { return x != x; }

int write_checksum(char * filename, DML_Checksum * checksum) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=1;
  char message[500];
  n_uint64_t bytes;
  /*   char * message; */


  if(g_cart_id == 0) {
    ofs = fopen(filename, "a");

    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }

    sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum><version>1.0</version><suma>%#010x</suma><sumb>%#010x</sumb></scidacChecksum>", (*checksum).suma, (*checksum).sumb);
    bytes = strlen( message );
    limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-checksum", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);

    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }

  return(0);
}

int write_inverter_info(const double epssq, const int iter, const int heavy, 
			const int append, char * filename, const int mms) {
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=1, MB_flag=1, status=0;
#ifdef MPI
  MPI_Status mpi_status;
#endif
  char message[5000];
  n_uint64_t bytes;
  struct timeval t1;
  

  gettimeofday(&t1,NULL);
  if(mms > -1) {
    sprintf(message,"\n result is for Q^dagger Q!\n multiple mass solver\n epssq = %e\n noiter = %d\n kappa = %f, inverted mu = %f, lowest mu = %f\n time = %ld\n hmcversion = %s\n date = %s", 
	    epssq, iter, g_kappa, g_extra_masses[mms]/2./g_kappa,
	    g_mu/2./g_kappa,t1.tv_sec, PACKAGE_VERSION, 
	    ctime(&t1.tv_sec));
  }
  else if(!heavy) {
    sprintf(message,"\n epssq = %e\n noiter = %d\n kappa = %f, mu = %f\n time = %ld\n hmcversion = %s\n date = %s", 
	    epssq, iter, g_kappa, g_mu/2./g_kappa,t1.tv_sec, PACKAGE_VERSION, 
	    ctime(&t1.tv_sec));
  }
  else {
    sprintf(message,"\n epssq = %e\n noiter = %d\n kappa = %f, mubar = %f, epsbar=%f\n time = %ld\n hmcversion = %s\n date = %s", 
	    epssq, iter, g_kappa, g_mubar/2./g_kappa, g_epsbar/2./g_kappa, t1.tv_sec, PACKAGE_VERSION, 
	    ctime(&t1.tv_sec));
  }
  bytes = strlen( message );
  if(g_cart_id == 0) {
    if(append) {
      ofs = fopen(filename, "a");
    }
    else ofs = fopen(filename, "w");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    
    limeheader = limeCreateHeader(MB_flag, ME_flag, "inverter-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);
    
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
    
  }
  return(0);
}

int write_xlf_info(const double plaq, const int counter, char * filename, 
		   const int append, char * data_buf) {
  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=1, MB_flag=1, status=0;
#ifdef MPI
  MPI_Status mpi_status;
#endif
  char message[5000];
  char * buf;
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
  if(data_buf != (char*)NULL) {
    bytes = strlen( data_buf );
    buf = data_buf;
  }
  else {
    bytes = strlen( message );
    buf = message;
  }
  if(g_cart_id == 0) {
    if(append) {
      ofs = fopen(filename, "a");
    }
    else ofs = fopen(filename, "w");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    
    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(buf, &bytes, limewriter);
    
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
    
  }
  return(0);
}

char * read_message(char * filename, char * type) {
  
  FILE * ifs;
  char * data_buf = NULL;
  int status;
  n_uint64_t bytes, read_bytes;
  char * header_type;
  LimeReader * limereader;

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n return empty message... (read_xlf_info)\n", filename);
    return(NULL);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n return empty message... (read_xlf_info)\n");
    return(NULL);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp(type, header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no %s record found in file %s\n", type, filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(NULL);
  }
  bytes = limeReaderBytes(limereader);

  data_buf = (char *)malloc(bytes+1);
  if( (data_buf = (char *)malloc(bytes+1)) == (char *)NULL) {
    fprintf(stderr, "Couldn't malloc data buf\n");
    return(NULL);
  }
  read_bytes = bytes;
  status = limeReaderReadData((void *)data_buf, &read_bytes, limereader);

  if( status < 0 ) {
    if( status != LIME_EOR ) {
      fprintf(stderr, "LIME read error occurred: status= %d  %llu bytes wanted, %llu read\n",
	      status, (unsigned long long)bytes,
	      (unsigned long long)read_bytes);
      free(data_buf);
      data_buf = NULL;
      return(NULL);
    }
  }
  data_buf[bytes]='\0';

  limeDestroyReader(limereader);
  fclose(ifs);

  return(data_buf);
}

int write_message(char * filename, char * data_buf, char * type, const int append) {

  FILE * ofs;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  /* Message end and Message begin flag */
  int ME_flag=1, MB_flag=1, status=0;
#ifdef MPI
  MPI_Status mpi_status;
#endif
  n_uint64_t bytes;

  if(data_buf == (char*)NULL) return(0);

  bytes = strlen( data_buf );

  if(g_cart_id == 0) {
    if(append) {
      ofs = fopen(filename, "a");
    }
    else ofs = fopen(filename, "w");
    if(ofs == (FILE*)NULL) {
      fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    
    limeheader = limeCreateHeader(MB_flag, ME_flag, type, bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header error %d\n", status);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(data_buf, &bytes, limewriter);
    
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
    
  }
  return(0);
}

int big_endian(){
  union{
    int l;
    char c[sizeof(int)];
  } u;

  u.l=1;
  return(u.c[sizeof(int) - 1] == 1);
}

void write_su3(su3 * up, FILE * f) {
  fprintf(f,"%f %f %f %f %f %f \n%f %f %f %f %f %f \n%f %f %f %f %f %f %d\n\n",
	     (*up).c00.re, (*up).c00.im, (*up).c01.re, (*up).c01.im,
	     (*up).c02.re, (*up).c02.im, (*up).c10.re, (*up).c10.im,
	     (*up).c11.re, (*up).c11.im, (*up).c12.re, (*up).c12.im,
	     (*up).c20.re, (*up).c20.im, (*up).c21.re, (*up).c21.im,
	     (*up).c22.re, (*up).c22.im, g_cart_id);
}

void byte_swap(void * ptr, int nmemb){
  int j;
  char char_in[4];
  char * in_ptr;
  int * int_ptr;

  for(j = 0, int_ptr = (int *) ptr; j < nmemb; j++, int_ptr++){
    in_ptr = (char *) int_ptr;
    
    char_in[0] = in_ptr[0];
    char_in[1] = in_ptr[1];
    char_in[2] = in_ptr[2];
    char_in[3] = in_ptr[3];

    in_ptr[0] = char_in[3];
    in_ptr[1] = char_in[2];
    in_ptr[2] = char_in[1];
    in_ptr[3] = char_in[0];
  }
}

void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr, * double_out_ptr;

  double_in_ptr = (double *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) double_in_ptr;
    char_out_ptr = (char *) double_out_ptr;
    
    char_out_ptr[7] = char_in_ptr[0];
    char_out_ptr[6] = char_in_ptr[1];
    char_out_ptr[5] = char_in_ptr[2];
    char_out_ptr[4] = char_in_ptr[3];
    char_out_ptr[3] = char_in_ptr[4];
    char_out_ptr[2] = char_in_ptr[5];
    char_out_ptr[1] = char_in_ptr[6];
    char_out_ptr[0] = char_in_ptr[7];
    double_in_ptr++;
    double_out_ptr++;
  }
}

void single2double(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) in_ptr;
  double * double_ptr = (double*) out_ptr;

  for(i = 0; i < nmemb; i++) {
    (*double_ptr) = (double) (*float_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void double2single(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) out_ptr;
  double * double_ptr = (double*) in_ptr;

  for(i = 0; i < nmemb; i++) {
    (*float_ptr) = (float) (*double_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_out_ptr;
  float * float_in_ptr;
  float tmp;

  float_in_ptr = (float *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  char_out_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) float_in_ptr;
    
    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];
    (*double_out_ptr) = (double) tmp;
    float_in_ptr++;
    double_out_ptr++;
  }
}

void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr;
  float * float_out_ptr;
  float tmp;

  float_out_ptr = (float *) out_ptr;
  double_in_ptr = (double *) in_ptr;
  char_in_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    tmp = (float) (*double_in_ptr);
    char_out_ptr = (char*) float_out_ptr;

    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];

    float_out_ptr++;
    double_in_ptr++;
  }
}


void single2double_cm(spinor * const R, float * const S) {
  (*R).s0.c0.re = (double) S[0];
  (*R).s0.c0.im = (double) S[1];
  (*R).s0.c1.re = (double) S[2];
  (*R).s0.c1.im = (double) S[3];
  (*R).s0.c2.re = (double) S[4];
  (*R).s0.c2.im = (double) S[5];
  (*R).s1.c0.re = (double) S[6];
  (*R).s1.c0.im = (double) S[7];
  (*R).s1.c1.re = (double) S[8];
  (*R).s1.c1.im = (double) S[9];
  (*R).s1.c2.re = (double) S[10];
  (*R).s1.c2.im = (double) S[11];
  (*R).s2.c0.re = (double) S[12];
  (*R).s2.c0.im = (double) S[13];
  (*R).s2.c1.re = (double) S[14];
  (*R).s2.c1.im = (double) S[15];
  (*R).s2.c2.re = (double) S[16];
  (*R).s2.c2.im = (double) S[17];
  (*R).s3.c0.re = (double) S[18];
  (*R).s3.c0.im = (double) S[19];
  (*R).s3.c1.re = (double) S[20];
  (*R).s3.c1.im = (double) S[21];
  (*R).s3.c2.re = (double) S[22];
  (*R).s3.c2.im = (double) S[23];
}

void double2single_cm(float * const S, spinor * const R) {
  S[0]  = (float) (*R).s0.c0.re ;
  S[1]  = (float) (*R).s0.c0.im ;
  S[2]  = (float) (*R).s0.c1.re ;
  S[3]  = (float) (*R).s0.c1.im ;
  S[4]  = (float) (*R).s0.c2.re ;
  S[5]  = (float) (*R).s0.c2.im ;
  S[6]  = (float) (*R).s1.c0.re ;
  S[7]  = (float) (*R).s1.c0.im ;
  S[8]  = (float) (*R).s1.c1.re ;
  S[9]  = (float) (*R).s1.c1.im ;
  S[10] = (float) (*R).s1.c2.re ;
  S[11] = (float) (*R).s1.c2.im ;
  S[12] = (float) (*R).s2.c0.re ;
  S[13] = (float) (*R).s2.c0.im ;
  S[14] = (float) (*R).s2.c1.re ;
  S[15] = (float) (*R).s2.c1.im ;
  S[16] = (float) (*R).s2.c2.re ;
  S[17] = (float) (*R).s2.c2.im ;
  S[18] = (float) (*R).s3.c0.re ;
  S[19] = (float) (*R).s3.c0.im ;
  S[20] = (float) (*R).s3.c1.re ;
  S[21] = (float) (*R).s3.c1.im ;
  S[22] = (float) (*R).s3.c2.re ;
  S[23] = (float) (*R).s3.c2.im ;
}

void zero_spinor(spinor * const R) {
  (*R).s0.c0.re = 0.;
  (*R).s0.c0.im = 0.;
  (*R).s0.c1.re = 0.;
  (*R).s0.c1.im = 0.;
  (*R).s0.c2.re = 0.;
  (*R).s0.c2.im = 0.;
  (*R).s1.c0.re = 0.;
  (*R).s1.c0.im = 0.;
  (*R).s1.c1.re = 0.;
  (*R).s1.c1.im = 0.;
  (*R).s1.c2.re = 0.;
  (*R).s1.c2.im = 0.;
  (*R).s2.c0.re = 0.;
  (*R).s2.c0.im = 0.;
  (*R).s2.c1.re = 0.;
  (*R).s2.c1.im = 0.;
  (*R).s2.c2.re = 0.;
  (*R).s2.c2.im = 0.;
  (*R).s3.c0.re = 0.;
  (*R).s3.c0.im = 0.;
  (*R).s3.c1.re = 0.;
  (*R).s3.c1.im = 0.;
  (*R).s3.c2.re = 0.;
  (*R).s3.c2.im = 0.;
}
