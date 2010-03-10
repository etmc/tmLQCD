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

#include"utils.ih"

int isnan_f  (float       x) { return x != x; }
int isnan_d  (double      x) { return x != x; }
int isnan_ld (long double x) { return x != x; }


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

#if BYTE_ORDER == LITTLE_ENDIAN

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

#else

void byte_swap(void * ptr, int nmemb){
  return;
}
#endif

#if BYTE_ORDER == LITTLE_ENDIAN

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
  return;
}

#else

void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb){
  memcpy(out_ptr, in_ptr, 8*nmemb);
  return;
}

#endif

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

#if BYTE_ORDER == LITTLE_ENDIAN

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
  return;
}

#else 

void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
  single2double(out_ptr, in_ptr, nmemb);
  return;
}

#endif

#if BYTE_ORDER == LITTLE_ENDIAN

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
  return;
}

#else

void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
  double2single(out_ptr, in_ptr, nmemb);
  return;
}
#endif

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
