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

#if (defined BGQ && defined XLC)
// r = w*r
// w a 6x6 colour matrix, r upper or lower half of a spinor
# define _mul_colourmatrix(w, r)		\
  r0 = vec_ld(0L, (double*) r);			\
  r1 = vec_ld(32L, (double*) r);		\
  r2 = vec_ld(64L, (double*) r);		\
  _Pragma("unroll(6)")				\
  for(int c0 = 0; c0 < 6; c0++) {		\
    w0 = vec_ld(0L, (double*) &w[6*c0]);	\
    w1 = vec_ld(32L, (double*) &w[6*c0]);	\
    w2 = vec_ld(64L, (double*) &w[6*c0]);	\
    r3 = vec_xmul(r0, w0);			\
    r4 = vec_xmul(r1, w1);			\
    r5 = vec_xmul(r2, w2);			\
    r3 = vec_xxnpmadd(w0, r0, r3);		\
    r4 = vec_xxnpmadd(w1, r1, r4);		\
    r5 = vec_xxnpmadd(w2, r2, r5);		\
    r6 = vec_add(r3, r4);			\
    r3 = vec_add(r5, r6);			\
    /* there must a smarter solution for this */	\
    r[c0] = r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
  }

#else
// needs temporary space _Complex double s[6]
# define _mul_colourmatrix(w, r)		\
  for(int c0 = 0; c0 < 6; c0++) {		\
    s[c0] = w[c0*6] * r[0];			\
    s[c0] += w[c0*6 + 1] * r[1];		\
    s[c0] += w[c0*6 + 2] * r[2];		\
    s[c0] += w[c0*6 + 3] * r[3];		\
    s[c0] += w[c0*6 + 4] * r[4];		\
    s[c0] += w[c0*6 + 5] * r[5];		\
  }						\
  memcpy(r, s, 6*sizeof(_Complex double));
#endif

#if (defined BGQ && defined XLC)
// r = w*r
// w a 6x6 colour matrix, r_s and r_c upper or lower halfs of a spinor
# define _mul_colourmatrix_nd(w, r_s, r_c)		\
  r0 = vec_ld(0L, (double*) r_s);			\
  r6 = vec_ld(0L, (double*) r_c);			\
  r1 = vec_ld(32L, (double*) r_s);			\
  r2 = vec_ld(64L, (double*) r_s);			\
  r7 = vec_ld(32L, (double*) r_c);			\
  r8 = vec_ld(64L, (double*) r_c);			\
  _Pragma("unroll(6)")					\
  for(int c0 = 0; c0 < 6; c0++) {			\
    w0 = vec_ld(0L, (double*) &w[6*c0]);		\
    w1 = vec_ld(32L, (double*) &w[6*c0]);		\
    w2 = vec_ld(64L, (double*) &w[6*c0]);		\
    r3  = vec_xmul(r0, w0);				\
    r4  = vec_xmul(r1, w1);				\
    r5  = vec_xmul(r2, w2);				\
    r9  = vec_xmul(r6, w0);				\
    r10 = vec_xmul(r7, w1);				\
    r11 = vec_xmul(r8, w2);				\
    r3  = vec_xxnpmadd(w0, r0, r3);			\
    r4  = vec_xxnpmadd(w1, r1, r4);			\
    r5  = vec_xxnpmadd(w2, r2, r5);			\
    r9  = vec_xxnpmadd(w0, r6, r9);			\
    r10 = vec_xxnpmadd(w1, r7, r10);			\
    r11 = vec_xxnpmadd(w2, r8, r11);			\
    r3 = vec_add(r3, r4);				\
    r9 = vec_add(r9, r10);				\
    r3 = vec_add(r5, r3);				\
    r9 = vec_add(r11, r9);				\
    /* there must a smarter solution for this */	\
    r_s[c0] = r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
    r_c[c0] = r9[0] + I*r9[1] + r9[2] + I*r9[3];	\
  }

#else
// needs temporary space _Complex double s[6]
# define _mul_colourmatrix_nd(w, r_s, r_c)		\
  for(int c0 = 0; c0 < 6; c0++) {			\
    s[c0] = w[c0*6] * r_s[0];				\
    s[c0+6] = w[c0*6] * r_c[0];				\
    s[c0] += w[c0*6 + 1] * r_s[1];			\
    s[c0+6] += w[c0*6 + 1] * r_c[1];			\
    s[c0] += w[c0*6 + 2] * r_s[2];			\
    s[c0+6] += w[c0*6 + 2] * r_c[2];			\
    s[c0] += w[c0*6 + 3] * r_s[3];			\
    s[c0+6] += w[c0*6 + 3] * r_c[3];			\
    s[c0] += w[c0*6 + 4] * r_s[4];			\
    s[c0+6] += w[c0*6 + 4] * r_c[4];			\
    s[c0] += w[c0*6 + 5] * r_s[5];			\
    s[c0+6] += w[c0*6 + 5] * r_c[5];			\
  }							\
  memcpy(r_s, s, 6*sizeof(_Complex double));		\
  memcpy(r_c, s+6, 6*sizeof(_Complex double));
#endif


#if (defined BGQ && defined XLC)
// s = w*r
// w a 6x6 colour matrix, r,s upper or lower half of a spinor
# define _mul_colourmatrix_assign(s, w, r)		\
  r0 = vec_ld(0L, (double*) r);			\
  r1 = vec_ld(32L, (double*) r);		\
  r2 = vec_ld(64L, (double*) r);		\
  _Pragma("unroll(6)")				\
  for(int c0 = 0; c0 < 6; c0++) {		\
    w0 = vec_ld(0L, (double*) &w[6*c0]);	\
    w1 = vec_ld(32L, (double*) &w[6*c0]);	\
    w2 = vec_ld(64L, (double*) &w[6*c0]);	\
    r3 = vec_xmul(r0, w0);			\
    r4 = vec_xmul(r1, w1);			\
    r5 = vec_xmul(r2, w2);			\
    r3 = vec_xxnpmadd(w0, r0, r3);		\
    r4 = vec_xxnpmadd(w1, r1, r4);		\
    r5 = vec_xxnpmadd(w2, r2, r5);		\
    r6 = vec_add(r3, r4);			\
    r3 = vec_add(r5, r6);			\
    s[c0] = r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
  }

#else

# define _mul_colourmatrix_assign(s, w, r)	\
  for(int c0 = 0; c0 < 6; c0++) {		\
    s[c0] = w[c0*6] * r[0];			\
    s[c0] += w[c0*6 + 1] * r[1];		\
    s[c0] += w[c0*6 + 2] * r[2];		\
    s[c0] += w[c0*6 + 3] * r[3];		\
    s[c0] += w[c0*6 + 4] * r[4];		\
    s[c0] += w[c0*6 + 5] * r[5];		\
  }

#endif

#if (defined BGQ && defined XLC)
// r = w*r
// w a 6x6 colour matrix, r_s and r_c upper or lower halfs of a spinor
# define _mul_colourmatrix_assign_nd(s_s, s_c, w, r_s, r_c)	\
  r0 = vec_ld(0L, (double*) r_s);			\
  r6 = vec_ld(0L, (double*) r_c);			\
  r1 = vec_ld(32L, (double*) r_s);			\
  r2 = vec_ld(64L, (double*) r_s);			\
  r7 = vec_ld(32L, (double*) r_c);			\
  r8 = vec_ld(64L, (double*) r_c);			\
  _Pragma("unroll(6)")					\
  for(int c0 = 0; c0 < 6; c0++) {			\
    w0 = vec_ld(0L, (double*) &w[6*c0]);		\
    w1 = vec_ld(32L, (double*) &w[6*c0]);		\
    w2 = vec_ld(64L, (double*) &w[6*c0]);		\
    r3  = vec_xmul(r0, w0);				\
    r4  = vec_xmul(r1, w1);				\
    r5  = vec_xmul(r2, w2);				\
    r9  = vec_xmul(r6, w0);				\
    r10 = vec_xmul(r7, w1);				\
    r11 = vec_xmul(r8, w2);				\
    r3  = vec_xxnpmadd(w0, r0, r3);			\
    r4  = vec_xxnpmadd(w1, r1, r4);			\
    r5  = vec_xxnpmadd(w2, r2, r5);			\
    r9  = vec_xxnpmadd(w0, r6, r9);			\
    r10 = vec_xxnpmadd(w1, r7, r10);			\
    r11 = vec_xxnpmadd(w2, r8, r11);			\
    r3 = vec_add(r3, r4);				\
    r9 = vec_add(r9, r10);				\
    r3 = vec_add(r5, r3);				\
    r9 = vec_add(r11, r9);				\
    /* there must a smarter solution for this */	\
    s_s[c0] = r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
    s_c[c0] = r9[0] + I*r9[1] + r9[2] + I*r9[3];	\
  }

#else

# define _mul_colourmatrix_assign_nd(s_s, s_c, w, r_s, r_c)	\
  for(int c0 = 0; c0 < 6; c0++) {				\
    s_s[c0] = w[c0*6] * r_s[0];					\
    s_c[c0] = w[c0*6] * r_c[0];					\
    s_s[c0] += w[c0*6 + 1] * r_s[1];				\
    s_c[c0] += w[c0*6 + 1] * r_c[1];				\
    s_s[c0] += w[c0*6 + 2] * r_s[2];				\
    s_c[c0] += w[c0*6 + 2] * r_c[2];				\
    s_s[c0] += w[c0*6 + 3] * r_s[3];				\
    s_c[c0] += w[c0*6 + 3] * r_c[3];				\
    s_s[c0] += w[c0*6 + 4] * r_s[4];				\
    s_c[c0] += w[c0*6 + 4] * r_c[4];				\
    s_s[c0] += w[c0*6 + 5] * r_s[5];				\
    s_c[c0] += w[c0*6 + 5] * r_c[5];				\
  }
#endif

#define _mul_colourmatrix_add_i_mul_assign(s, w, mu, r)		\
  for(int c0 = 0; c0 < 6; c0++) {				\
    s[c0] = w[c0*6] * r[0];					\
    s[c0] += w[c0*6 + 1] * r[1];				\
    s[c0] += w[c0*6 + 2] * r[2];				\
    s[c0] += w[c0*6 + 3] * r[3];				\
    s[c0] += w[c0*6 + 4] * r[4];				\
    s[c0] += w[c0*6 + 5] * r[5];				\
    s[c0] += I*mu*r[c0];				\
  }

#if (defined BGQ && defined XLC)
#define _mul_colourmatrix_add_i_mul_sub_assign(s, w, mu, r, t)	\
  r0 = vec_ld(0L, (double*) r);			\
  r1 = vec_ld(32L, (double*) r);		\
  r2 = vec_ld(64L, (double*) r);		\
  r6 = vec_ld(0L, (double*) t);			\
  r7 = vec_ld(32L, (double*) t);		\
  r8 = vec_ld(64L, (double*) t);		\
  tmp = vec_splats(mu);				\
  r6 = vec_neg(r6);				\
  r7 = vec_neg(r7);				\
  r8 = vec_neg(r8);				\
  r3 = vec_xxnpmadd(r0, tmp, r6);		\
  r4 = vec_xxnpmadd(r1, tmp, r7);		\
  r5 = vec_xxnpmadd(r2, tmp, r8);		\
  vec_st(r3, 0L, (double*) s);			\
  vec_st(r4, 32L, (double*) s);			\
  vec_st(r5, 64L, (double*) s);			\
  _Pragma("unroll(6)")				\
  for(int c0 = 0; c0 < 6; c0++) {		\
    w0 = vec_ld(0L, (double*) &w[6*c0]);	\
    w1 = vec_ld(32L, (double*) &w[6*c0]);	\
    w2 = vec_ld(64L, (double*) &w[6*c0]);	\
    r3 = vec_xmul(r0, w0);			\
    r4 = vec_xmul(r1, w1);			\
    r5 = vec_xmul(r2, w2);			\
    r3 = vec_xxnpmadd(w0, r0, r3);		\
    r4 = vec_xxnpmadd(w1, r1, r4);		\
    r5 = vec_xxnpmadd(w2, r2, r5);		\
    r3 = vec_add(r3, r4);			\
    r3 = vec_add(r5, r3);			\
    s[c0] += r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
  }

#else

#define _mul_colourmatrix_add_i_mul_sub_assign(s, w, mu, r, t)	\
  for(int c0 = 0; c0 < 6; c0++) {				\
    s[c0] = w[c0*6] * r[0];					\
    s[c0] += w[c0*6 + 1] * r[1];				\
    s[c0] += w[c0*6 + 2] * r[2];				\
    s[c0] += w[c0*6 + 3] * r[3];				\
    s[c0] += w[c0*6 + 4] * r[4];				\
    s[c0] += w[c0*6 + 5] * r[5];				\
    s[c0] += I*mu*r[c0] - t[c0];				\
  }
#endif

#if (defined BGQ && defined XLC)
#define _mul_colourmatrix_add_i_mul_sub_nassign(s, w, mu, r, t)	\
  r0 = vec_ld(0L, (double*) r);			\
  r1 = vec_ld(32L, (double*) r);		\
  r2 = vec_ld(64L, (double*) r);		\
  r6 = vec_ld(0L, (double*) t);			\
  r7 = vec_ld(32L, (double*) t);		\
  r8 = vec_ld(64L, (double*) t);		\
  tmp = vec_splats(mu);				\
  r0 = vec_neg(r0);				\
  r1 = vec_neg(r1);				\
  r2 = vec_neg(r2);				\
  r3 = vec_xxnpmadd(r0, tmp, r6);		\
  r4 = vec_xxnpmadd(r1, tmp, r7);		\
  r5 = vec_xxnpmadd(r2, tmp, r8);		\
  vec_st(r3, 0L, (double*) s);			\
  vec_st(r4, 32L, (double*) s);			\
  vec_st(r5, 64L, (double*) s);			\
  _Pragma("unroll(6)")				\
  for(int c0 = 0; c0 < 6; c0++) {		\
    w0 = vec_ld(0L, (double*) &w[6*c0]);	\
    w1 = vec_ld(32L, (double*) &w[6*c0]);	\
    w2 = vec_ld(64L, (double*) &w[6*c0]);	\
    r3 = vec_xmul(r0, w0);			\
    r4 = vec_xmul(r1, w1);			\
    r5 = vec_xmul(r2, w2);			\
    r3 = vec_xxnpmadd(w0, r0, r3);		\
    r4 = vec_xxnpmadd(w1, r1, r4);		\
    r5 = vec_xxnpmadd(w2, r2, r5);		\
    r3 = vec_add(r3, r4);			\
    r3 = vec_add(r5, r3);			\
    s[c0] += r3[0] + I*r3[1] + r3[2] + I*r3[3];	\
  }

#else

#define _mul_colourmatrix_add_i_mul_sub_nassign(s, w, mu, r, t)	\
  for(int c0 = 0; c0 < 6; c0++) {				\
    s[c0] = -w[c0*6] * r[0];					\
    s[c0] -= w[c0*6 + 1] * r[1];				\
    s[c0] -= w[c0*6 + 2] * r[2];				\
    s[c0] -= w[c0*6 + 3] * r[3];				\
    s[c0] -= w[c0*6 + 4] * r[4];				\
    s[c0] -= w[c0*6 + 5] * r[5];				\
    s[c0] -= I*mu*r[c0] - t[c0];				\
  }
#endif

#define _add_nd_massterm_sub(r_s, r_c, mubar, epsbar, s_s, s_c, t_s, t_c) \
  for(int c0 = 0; c0 < 6; c0++) {					\
    r_s[c0] += +I*mubar*s_s[c0] + epsbar*s_c[c0] - t_s[c0];		\
    r_c[c0] += -I*mubar*s_c[c0] + epsbar*s_s[c0] - t_c[c0];		\
  }

#define _add_nd_massterm(r_s, r_c, mubar, epsbar, s_s, s_c)	\
  for(int c0 = 0; c0 < 6; c0++) {				\
    r_s[c0] += +I*mubar*s_s[c0] + epsbar*s_c[c0];		\
    r_c[c0] += -I*mubar*s_c[c0] + epsbar*s_s[c0];		\
  }

#define _add_nd_massterm_nsub(r_s, r_c, mubar, epsbar, s_s, s_c, t_s, t_c) \
  for(int c0 = 0; c0 < 6; c0++) {					\
    r_s[c0] = -(r_s[c0] + I*mubar*s_s[c0] + epsbar*s_c[c0] - t_s[c0]);	\
    r_c[c0] = -(r_c[c0] - I*mubar*s_c[c0] + epsbar*s_s[c0] - t_c[c0]);	\
  }
