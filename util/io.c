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

/****************************************************
 * IO routines:
 *
 * read_lime_gauge_field_doubleprec
 *
 * read_lime_gauge_field_singleprec
 *
 * Autor: 
 *        Carsten Urbach <urbach@ifh.de>
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
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<sys/time.h> 
#include<sys/types.h>
#include "io.h"

#define MAXBUF 1048576

void byte_swap(void *ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_singleprec(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb);
void single2double(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb);
void double2single(void * out_ptr, void * in_ptr, int nmemb);
int big_endian();

int read_lime_gauge_field_doubleprec(double * config, char * filename,
				     const int T, const int LX, const int LY, const int LZ) {
  FILE * ifs;
  int t, x, y, z, status, p=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  double tmp[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*LY*LZ*T*72*sizeof(double)) {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
    exit(501);
  }

  bytes = (n_uint64_t)72*sizeof(double);

  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
	for(x = 0; x < LX; x++) {
	  p = (((t*LZ+z)*LY+y)*LX+x)*72;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(&config[p], tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(&config[p], &bytes, limereader);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int read_lime_gauge_field_singleprec(float * config, char * filename,
				     const int T, const int LX, const int LY, const int LZ){
  FILE * ifs;
  int t, x, y, z, status, p=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  float tmp[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*LY*LZ*T*72*sizeof(float)) {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
    exit(501);
  }

  bytes = (n_uint64_t)72*sizeof(float);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
	for(x = 0; x < LX; x++) {
	  p = (((t*LZ+z)*LY+y)*LX+x)*72;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign_singleprec(&config[p], tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(&config[p], &bytes, limereader);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
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

void byte_swap(void * ptr, int nmemb){
  int j;
  char char_in[4];
  char * in_ptr;
  int * int_ptr;

  for(j = 0, int_ptr = (int *) ptr; j < nmemb; j++, int_ptr++) {
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

void byte_swap_assign_singleprec(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  float * float_in_ptr, * float_out_ptr;

  float_in_ptr = (float *) in_ptr;
  float_out_ptr = (float *) out_ptr;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) float_in_ptr;
    char_out_ptr = (char *) float_out_ptr;
    
    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];
    float_in_ptr++;
    float_out_ptr++;
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

