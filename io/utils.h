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

#ifndef _UTILS_H
#define _UTILS_H

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <endian.h>
#include <string.h>

#include "su3.h"
#include <io/selector.h>
#include <io/params.h>
#include <io/dml.h>


/* These are factory functions, since the constructors for c-lime and lemon are different
   and they need different ways of opening files. Moving this to utility functions unclutters
   the main code, since we don't need additional #ifdefs anymore.
   Since lemon is collective and c-lime isn't, some care needs to be taken. For now, a writer
   will only be constructed on the node with g_cart_id == 0 for c-lime, while a reader will
   be created everywhere to exploit trivial parallellization. Be careful not to call
   construct_writer and friends from within a "if (node_num == 0)" type statement, because
   it will cause lemon to get deadlocked! */
void construct_writer(WRITER ** writer, char * filename, const int append);
void destruct_writer(WRITER * writer);

void construct_reader(READER ** reader, char * filename);
void destruct_reader(READER * reader);

void kill_with_error(LIME_FILE *fh, int const rank, char const *error);

int read_message(READER *reader, char **buffer);
int write_message(WRITER * writer, char const *buffer, uint64_t bytes);
void write_header(WRITER * writer, int MB, int ME, char const *type, uint64_t bytes);

void write_checksum(WRITER *writer, DML_Checksum const *checksum, char const *name);
void write_xlf_info(WRITER *writer, paramsXlfInfo const *info);
void write_xlf_info_xml(WRITER *writer, paramsXlfInfo const *info);
void write_inverter_info(WRITER * writer, paramsInverterInfo const *info);

void close_reader_record(READER *reader);
void close_writer_record(WRITER *writer);

void engineering(char *result, double value, char const *units);
int parse_checksum_xml(char *message, DML_Checksum *checksum);

int big_endian();
int write_ildg_format_xml(char *filename, LimeWriter * limewriter, const int precision);
void single2double_cm(spinor * const R, float * const S);
void double2single_cm(float * const S, spinor * const R);
void zero_spinor(spinor * const R);

int write_first_messages(FILE * parameterfile, char const * const executable, char const * const git_hash);
int parse_propagator_type(READER * reader);

int parse_ildgformat_xml(char *message, paramsIldgFormat *ildgformat);

inline static void byte_swap(void * ptr, int nmemb){
  int j;
  char char_in[8];
  char * in_ptr;
  double * d_ptr;

  for(j = 0, d_ptr = (double *) ptr; j < nmemb; j++, d_ptr++){
    in_ptr = (char *) d_ptr;
    
    char_in[0] = in_ptr[0];
    char_in[1] = in_ptr[1];
    char_in[2] = in_ptr[2];
    char_in[3] = in_ptr[3];
    char_in[4] = in_ptr[4];
    char_in[5] = in_ptr[5];
    char_in[6] = in_ptr[6];
    char_in[7] = in_ptr[7];

    in_ptr[0] = char_in[7];
    in_ptr[1] = char_in[6];
    in_ptr[2] = char_in[5];
    in_ptr[3] = char_in[4];
    in_ptr[4] = char_in[3];
    in_ptr[5] = char_in[2];
    in_ptr[6] = char_in[1];
    in_ptr[7] = char_in[0];
  }
}

inline static void byte_swap32(void * ptr, int nmemb){
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

inline static void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb){
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
  return;
}

inline static void byte_swap_assign32(void * out_ptr, void * in_ptr, int nmemb){
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
  return;
}


#if BYTE_ORDER == LITTLE_ENDIAN

inline static void be_to_cpu_assign(void * out_ptr, void * in_ptr, int nmemb){
  byte_swap_assign(out_ptr, in_ptr, nmemb);
  return;
}

#else

inline static void be_to_cpu_assign(void * out_ptr, void * in_ptr, int nmemb){
  memcpy(out_ptr, in_ptr, 8*nmemb);
  return;
}

#endif

inline static void single2double(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) in_ptr;
  double * double_ptr = (double*) out_ptr;

  for(i = 0; i < nmemb; i++) {
    (*double_ptr) = (double) (*float_ptr);

    float_ptr++;
    double_ptr++;
  }

}

inline static void double2single(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) out_ptr;
  double * double_ptr = (double*) in_ptr;

  for(i = 0; i < nmemb; i++) {
    (*float_ptr) = (float) (*double_ptr);

    float_ptr++;
    double_ptr++;
  }

}

#if BYTE_ORDER == LITTLE_ENDIAN

inline static void be_to_cpu_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
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
  return;
}

#else 

inline static void be_to_cpu_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
  single2double(out_ptr, in_ptr, nmemb);
  return;
}

#endif

#if BYTE_ORDER == LITTLE_ENDIAN

inline static void be_to_cpu_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
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
  return;
}

#else

inline static void be_to_cpu_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
  double2single(out_ptr, in_ptr, nmemb);
  return;
}
#endif


#endif
