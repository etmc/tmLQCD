/***********************************************************************
 * Copyright (C) 2001 Martin Hasenbusch
 * Modified by Bartosz Kostrzewa (2012 void versions)
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
 *
 * File expo.c
 *
 *
 * The externally accessible functions are
 * 
 * void exposu3(su3* const vr, const su3adj* const p);
 * extern void exposu3_check(su3* const vr, const su3adj* const p, int im);
 * extern void restoresu3(su3* const vr,const su3* const u);
 * extern void restoresu3_in_place(su3* const u);
 * extern void exposu3_in_place(su3* const u);
 *
 * Author: Martin Hasenbusch <martin.hasenbusch@desy.de>
 * Tue Aug 28 10:06:56 MEST 2001
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef SSE
# undef SSE
#endif
#ifdef SSE2
# undef SSE2
#endif
#ifdef SSE3
# undef SSE3
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "expo.h"

void exposu3(su3* const vr, const su3adj* const p) {
  int i;
  su3 ALIGN v,v2;
  double ALIGN fac,r;
  double ALIGN a,b;
  _Complex double ALIGN a0,a1,a2,a1p;

  /* it writes 'p=vec(h_{j,mu})' in matrix form 'v' */  
  _make_su3(v,*p);
  /* calculates v^2 */
  _su3_times_su3(v2,v,v);
  /* */
  a = 0.5 * (creal(v2.c00) + creal(v2.c11) + creal(v2.c22));
  /* 1/3 imaginary part of tr v*v2 */
  b = 0.33333333333333333 * cimag(v.c00 * v2.c00 + v.c01 * v2.c10 + v.c02 * v2.c20 +
                                  v.c10 * v2.c01 + v.c11 * v2.c11 + v.c12 * v2.c21 +
                                  v.c20 * v2.c02 + v.c21 * v2.c12 + v.c22 * v2.c22  );
  a0  = 0.16059043836821615e-9;
  a1  = 0.11470745597729725e-10;
  a2  = 0.76471637318198165e-12;
  fac = 0.20876756987868099e-8;      /*  1/12! */
  r   = 12.0;
  for(i = 3; i <= 15; ++i)
  {
    a1p = a0 + a * a2;
    a0 = fac + b * I * a2;
    a2 = a1;
    a1 = a1p;
    fac *= r;
    r -= 1.0;
  }
  /* vr = a0 + a1*v + a2*v2 */
  vr->c00 = a0 + a1 * v.c00 + a2 * v2.c00;
  vr->c01 =      a1 * v.c01 + a2 * v2.c01;
  vr->c02 =      a1 * v.c02 + a2 * v2.c02;
  vr->c10 =      a1 * v.c10 + a2 * v2.c10;
  vr->c11 = a0 + a1 * v.c11 + a2 * v2.c11;
  vr->c12 =      a1 * v.c12 + a2 * v2.c12;
  vr->c20 =      a1 * v.c20 + a2 * v2.c20;
  vr->c21 =      a1 * v.c21 + a2 * v2.c21;
  vr->c22 = a0 + a1 * v.c22 + a2 * v2.c22;
}

void exposu3_check(su3* const vr, const su3adj* const p, int im) {
  /* compute the result by taylor series */
  su3 ALIGN v,v2,v3;
  double ALIGN fac;
  int i;

  _make_su3(v, *p);
  _su3_one(*vr);
  _su3_acc(*vr, v); 
  _su3_times_su3(v2, v, v);
  _su3_refac_acc(*vr, 0.5, v2);
  fac = 0.5;
  for(i = 3; i <= im; i++) {
    fac = fac/i;
    _su3_times_su3(v3, v2, v);
    _su3_refac_acc(*vr, fac, v3); 
    _su3_assign(v2, v3); 
  }
}

void restoresu3(su3* const vr, const su3* const u) {
  double ALIGN n0,n1;

  /* normalize rows 1 and 2 */
  n0 = 1.0 / sqrt(conj(u->c00) * u->c00 + conj(u->c01) * u->c01 + conj(u->c02) * u->c02);
  n1 = 1.0 / sqrt(conj(u->c10) * u->c10 + conj(u->c11) * u->c11 + conj(u->c12) * u->c12);

  vr->c00 = n0 * u->c00;
  vr->c01 = n0 * u->c01;
  vr->c02 = n0 * u->c02;

  vr->c10 = n1 * u->c10;
  vr->c11 = n1 * u->c11;
  vr->c12 = n1 * u->c12;

  /* compute  row 3 as the conjugate of the cross-product of 1 and 2 */ 
  vr->c20 = conj(vr->c01 * vr->c12 - vr->c02 * vr->c11);
  vr->c21 = conj(vr->c02 * vr->c10 - vr->c00 * vr->c12);
  vr->c22 = conj(vr->c00 * vr->c11 - vr->c01 * vr->c10);
}

void restoresu3_in_place(su3* const u) {
  double ALIGN n0,n1;
    
  /* normalize rows 1 and 2 */
  n0 = 1.0 / sqrt(conj(u->c00) * u->c00 + conj(u->c01) * u->c01 + conj(u->c02) * u->c02);
  n1 = 1.0 / sqrt(conj(u->c10) * u->c10 + conj(u->c11) * u->c11 + conj(u->c12) * u->c12);
          
  u->c00 = n0 * u->c00;
  u->c01 = n0 * u->c01;
  u->c02 = n0 * u->c02;
                
  u->c10 = n1 * u->c10;
  u->c11 = n1 * u->c11;
  u->c12 = n1 * u->c12;
                      
  /* compute  row 3 as the conjugate of the cross-product of 1 and 2 */
  u->c20 = conj(u->c01 * u->c12 - u->c02 * u->c11);
  u->c21 = conj(u->c02 * u->c10 - u->c00 * u->c12);
  u->c22 = conj(u->c00 * u->c11 - u->c01 * u->c10);
}
                                
/* Exponentiates a hermitian 3x3 matrix Q */
/* Convenience function -- wrapper around Hasenbusch's implementation */
void exposu3_in_place(su3* const u) {
  su3adj ALIGN p;

  _trace_lambda(p, *u); /* Projects onto the Gell-Mann matrices */
  /* -2.0 to get su3 to su3adjoint consistency ****/
  p.d1 *= -0.5; p.d2 *= -0.5; p.d3 *= -0.5; p.d4 *= -0.5;
  p.d5 *= -0.5; p.d6 *= -0.5; p.d7 *= -0.5; p.d8 *= -0.5;
  exposu3(u,&p);
}
