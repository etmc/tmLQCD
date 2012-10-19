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

#ifndef _CLOVER_LEAF_H
#define _CLOVER_LEAF_H
#include "su3.h"
#include "hamiltonian_field.h"

extern su3 ** swm, ** swp;
extern const double tiny_t;

void sw_term(const su3 ** const gf, const double kappa, const double c_sw);
double sw_trace(const int ieo, const double mu);
double sw_trace_nd(const int ieo, const double mu, const double eps);
void sw_invert(const int ieo, const double mu);
void sw_invert_nd(const double mshift);
void sw_deriv(const int ieo, const double mu);
void sw_deriv_nd(const int ieo);
void sw_spinor(const int ieo, const spinor * const kk, const spinor * const ll, const double fac);
void sw_all(hamiltonian_field_t * const hf, const double kappa, const double c_sw);
int init_swpm(const int V);

void mult_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]);
void add_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]);
void sub_6x6(_Complex double a[6][6], _Complex double b[6][6], _Complex double d[6][6]);
void copy_6x6(_Complex double a[6][6], const _Complex double b[6][6]);

/*definitions needed for the functions sw_trace(int ieo) and sw_trace(int ieo)*/
inline void populate_6x6_matrix(_Complex double a[6][6], const su3 * const C, const int row, const int col) {
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

inline void get_3x3_block_matrix(su3 * const C, _Complex double a[6][6], const int row, const int col) {
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

// This function computes the trace-log part of the clover term
// in case of even/odd preconditioning
//
// it is expected that sw_term is called beforehand such that
// the array sw is populated properly

inline void add_tm(_Complex double a[6][6], const double mu) {
  for(int i = 0; i < 6; i++) {
    a[i][i] += I*mu;
  }
  return;
}

inline void add_shift_6x6(_Complex double a[6][6], const double mshift) {
  for(int i = 0; i < 6; i++) {
    a[i][i] += mshift;
  }
  return;
}


#endif
