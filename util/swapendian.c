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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>


void usage(){
  fprintf(stdout, "Usage:   swapendian [options]\n");
  fprintf(stdout, "Options: [-i input-filename]\n");
  fprintf(stdout, "         [-o output-filename]\n");
  fprintf(stdout, "         [-e single precision]\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}
void byte_swap_assign_singleprec(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);

int main(int argc,char *argv[]) {

  int c;
  FILE *ifs, *ofs;
  char * ifilename = NULL;
  char * ofilename = NULL;
  int single = 0;
  double tmpd, swapd;
  float tmps, swaps;
  int cnt = 0;

  while ((c = getopt(argc, argv, "h?i:o:e")) != -1) {
    switch (c) {
    case 'i': 
      ifilename = (char*)calloc(200, sizeof(char));
      strcpy(ifilename,optarg);
      break;
    case 'o':
      ofilename = (char*)calloc(200, sizeof(char));
      strcpy(ofilename,optarg);
      break;
    case 'e':
      single = 1;
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if(ifilename == NULL){
    fprintf(stderr, "input filename missing! Aborting...\n");
    exit(-1);
  }
  ifs = fopen(ifilename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", ifilename);
    exit(500);
  }

  if(ofilename == NULL){
    fprintf(stderr, "output filename missing! Aborting...\n");
    exit(-2);
  } 
  ofs = fopen(ofilename, "w");
  if(ofs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", ofilename);
    exit(500);
  }

  while(!feof(ifs)) {
    if(!single) {
      fread(&tmpd, sizeof(double), 1, ifs);
      if(!feof(ifs)) {
	cnt++;
	byte_swap_assign(&swapd, &tmpd, 1);
	fwrite(&swapd, sizeof(double), 1, ofs);
      }
    }
    else {
      fread(&tmps, sizeof(float), 1, ifs);
      if(!feof(ifs)) {
	cnt++;
	byte_swap_assign_singleprec(&swaps, &tmps, 1);
	fwrite(&swaps, sizeof(float), 1, ofs);
      }
    }
  }

  printf("Swapped endian for %d words\n", cnt);
  if(single) {
    printf("in- and output file in single precision\n");
  }
  else {
    printf("in- and output file in double precision\n");
  }

  fclose(ofs);
  fclose(ifs);

  return(0);
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
