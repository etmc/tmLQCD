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
	     creal(up->c00), cimag(up->c00), creal(up->c01), cimag(up->c01),
	     creal(up->c02), cimag(up->c02), creal(up->c10), cimag(up->c10),
	     creal(up->c11), cimag(up->c11), creal(up->c12), cimag(up->c12),
	     creal(up->c20), cimag(up->c20), creal(up->c21), cimag(up->c21),
	     creal(up->c22), cimag(up->c22), g_cart_id);
}



void single2double_cm(spinor * const R, float * const S) {
  R->s0.c0 = S[ 0] + S[ 1] * I;
  R->s0.c1 = S[ 2] + S[ 3] * I;
  R->s0.c2 = S[ 4] + S[ 5] * I;
  R->s1.c0 = S[ 6] + S[ 7] * I;
  R->s1.c1 = S[ 8] + S[ 9] * I;
  R->s1.c2 = S[10] + S[11] * I;
  R->s2.c0 = S[12] + S[13] * I;
  R->s2.c1 = S[14] + S[15] * I;
  R->s2.c2 = S[16] + S[17] * I;
  R->s3.c0 = S[18] + S[19] * I;
  R->s3.c1 = S[20] + S[21] * I;
  R->s3.c2 = S[22] + S[23] * I;
}

void double2single_cm(float * const S, spinor * const R) {
  S[ 0] = creal(R->s0.c0);
  S[ 1] = cimag(R->s0.c0);
  S[ 2] = creal(R->s0.c1);
  S[ 3] = cimag(R->s0.c1);
  S[ 4] = creal(R->s0.c2);
  S[ 5] = cimag(R->s0.c2);
  S[ 6] = creal(R->s1.c0);
  S[ 7] = cimag(R->s1.c0);
  S[ 8] = creal(R->s1.c1);
  S[ 9] = cimag(R->s1.c1);
  S[10] = creal(R->s1.c2);
  S[11] = cimag(R->s1.c2);
  S[12] = creal(R->s2.c0);
  S[13] = cimag(R->s2.c0);
  S[14] = creal(R->s2.c1);
  S[15] = cimag(R->s2.c1);
  S[16] = creal(R->s2.c2);
  S[17] = cimag(R->s2.c2);
  S[18] = creal(R->s3.c0);
  S[19] = cimag(R->s3.c0);
  S[20] = creal(R->s3.c1);
  S[21] = cimag(R->s3.c1);
  S[22] = creal(R->s3.c2);
  S[23] = cimag(R->s3.c2);
}

void zero_spinor(spinor * const R) {
  R->s0.c0 = 0.;
  R->s0.c1 = 0.;
  R->s0.c2 = 0.;
  R->s1.c0 = 0.;
  R->s1.c1 = 0.;
  R->s1.c2 = 0.;
  R->s2.c0 = 0.;
  R->s2.c1 = 0.;
  R->s2.c2 = 0.;
  R->s3.c0 = 0.;
  R->s3.c1 = 0.;
  R->s3.c2 = 0.;
}
