/***********************************************************************
 *
 * Copyright (C) 2005 Martin Hasenbusch
 *               2011 Carsten Urbach
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

/*definitions needed for the functions sw_trace(int ieo) and sw_trace_nd(int ieo)*/
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

static inline void six_mul_six(_Complex double c[6][6], _Complex double a[6][6], _Complex double b[6][6]) {
  for(unsigned int i = 0; i < 6; ++i) {
    for(unsigned int j = 0; j < 6; ++j) {
      c[i][j] = 0;
      for(unsigned int k = 0; k < 6; ++k) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return;
}

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

