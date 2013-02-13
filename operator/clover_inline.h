/***********************************************************************
 *
 * Copyright (C) 2011 Carsten Urbach
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

/*definitions needed for the functions sw_trace(int ieo) and sw_trace(int ieo)*/
static inline void populate_6x6_matrix(_Complex double a[6][6], const su3 * const C, const int row, const int col) {
  a[0+row][0+col] = C->c00;
  a[0+row][1+col] = C->c01;
  a[0+row][2+col] = C->c02;
  a[1+row][0+col] = C->c10;
  a[1+row][1+col] = C->c11;
  a[1+row][2+col] = C->c12;
  a[2+row][0+col] = C->c20;
  a[2+row][1+col] = C->c21;
  a[2+row][2+col] = C->c22;
  return;
}

static inline void populate_6x6_matrix2(_Complex double * a, const su3 * const C, const int row, const int col) {
  a[(0+row)*6 + 0+col] = C->c00;
  a[(0+row)*6 + 1+col] = C->c01;
  a[(0+row)*6 + 2+col] = C->c02;
  a[(1+row)*6 + 0+col] = C->c10;
  a[(1+row)*6 + 1+col] = C->c11;
  a[(1+row)*6 + 2+col] = C->c12;
  a[(2+row)*6 + 0+col] = C->c20;
  a[(2+row)*6 + 1+col] = C->c21;
  a[(2+row)*6 + 2+col] = C->c22;
  return;
}

static inline void populate_6x6_hc_matrix(_Complex double * a, const su3 * const C, const int row, const int col) {
  a[(0+row)*6 + 0+col] = conj(C->c00);
  a[(0+row)*6 + 1+col] = conj(C->c10);
  a[(0+row)*6 + 2+col] = conj(C->c20);
  a[(1+row)*6 + 0+col] = conj(C->c01);
  a[(1+row)*6 + 1+col] = conj(C->c11);
  a[(1+row)*6 + 2+col] = conj(C->c21);
  a[(2+row)*6 + 0+col] = conj(C->c02);
  a[(2+row)*6 + 1+col] = conj(C->c12);
  a[(2+row)*6 + 2+col] = conj(C->c22);
  return;
}


static inline void get_3x3_block_matrix(su3 * const C, _Complex double a[6][6], const int row, const int col) {
  C->c00 = a[0+row][0+col];
  C->c01 = a[0+row][1+col];
  C->c02 = a[0+row][2+col];
  C->c10 = a[1+row][0+col];
  C->c11 = a[1+row][1+col];
  C->c12 = a[1+row][2+col];
  C->c20 = a[2+row][0+col];
  C->c21 = a[2+row][1+col];
  C->c22 = a[2+row][2+col];
  return;
}

static inline void colour_plus_colour_addto_su3(su3 * const C, 
						_Complex double * a, _Complex double * b, 
						const double fac, const int row, const int col) {
  C->c00 += fac*( a[(0+row)*6 + 0+col] + b[(0+row)*6 + 0+col] );
  C->c01 += fac*( a[(0+row)*6 + 1+col] + b[(0+row)*6 + 1+col] );
  C->c02 += fac*( a[(0+row)*6 + 2+col] + b[(0+row)*6 + 2+col] );
  C->c10 += fac*( a[(1+row)*6 + 0+col] + b[(1+row)*6 + 0+col] );
  C->c11 += fac*( a[(1+row)*6 + 1+col] + b[(1+row)*6 + 1+col] );
  C->c12 += fac*( a[(1+row)*6 + 2+col] + b[(1+row)*6 + 2+col] );
  C->c20 += fac*( a[(2+row)*6 + 0+col] + b[(2+row)*6 + 0+col] );
  C->c21 += fac*( a[(2+row)*6 + 1+col] + b[(2+row)*6 + 1+col] );
  C->c22 += fac*( a[(2+row)*6 + 2+col] + b[(2+row)*6 + 2+col] );
}

static inline void colour_minus_colour_addto_su3(su3 * const C, 
						 _Complex double * a, _Complex double * b, 
						 const double fac, const int row, const int col) {
  C->c00 += fac*( a[(0+row)*6 + 0+col] - b[(0+row)*6 + 0+col] );
  C->c01 += fac*( a[(0+row)*6 + 1+col] - b[(0+row)*6 + 1+col] );
  C->c02 += fac*( a[(0+row)*6 + 2+col] - b[(0+row)*6 + 2+col] );
  C->c10 += fac*( a[(1+row)*6 + 0+col] - b[(1+row)*6 + 0+col] );
  C->c11 += fac*( a[(1+row)*6 + 1+col] - b[(1+row)*6 + 1+col] );
  C->c12 += fac*( a[(1+row)*6 + 2+col] - b[(1+row)*6 + 2+col] );
  C->c20 += fac*( a[(2+row)*6 + 0+col] - b[(2+row)*6 + 0+col] );
  C->c21 += fac*( a[(2+row)*6 + 1+col] - b[(2+row)*6 + 1+col] );
  C->c22 += fac*( a[(2+row)*6 + 2+col] - b[(2+row)*6 + 2+col] );
}

// This function computes the trace-log part of the clover term
// in case of even/odd preconditioning
//
// it is expected that sw_term is called beforehand such that
// the array sw is populated properly

static inline void add_tm(_Complex double a[6][6], const double mu) {
  for(int i = 0; i < 6; i++) {
    a[i][i] += I*mu;
  }
  return;
}

static inline void add_shift_6x6(_Complex double a[6][6], const double mshift) {
  for(int i = 0; i < 6; i++) {
    a[i][i] += mshift;
  }
  return;
}

