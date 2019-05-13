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
#include "float.h"
#include "global.h"

static double imag_det(const su3adj* p) {
  double d,tos3,o3,os3;
  tos3=2.0/sqrt(3.0);
  o3=1.0/3.0;
  os3=1.0/sqrt(3.0);
  
  d=tos3*(*p).d8*(o3*(*p).d8*(*p).d8-(*p).d3*(*p).d3)+2*((*p).d2*(*p).d4*(*p).d7-(*p).d1*(*p).d4*(*p).d6-(*p).d2*(*p).d5*(*p).d6-(*p).d1*(*p).d5*(*p).d7);
  d+=(os3*(*p).d8-(*p).d3)*((*p).d4*(*p).d4+(*p).d5*(*p).d5)+(os3*(*p).d8+(*p).d3)*((*p).d6*(*p).d6+(*p).d7*(*p).d7)-tos3*(*p).d8*((*p).d1*(*p).d1+(*p).d2*(*p).d2);	
  return d;
}

static void mul_su3alg(su3adj* p,double d) {
  (*p).d1*=d;
  (*p).d2*=d;
  (*p).d3*=d;
  (*p).d4*=d;
  (*p).d5*=d;
  (*p).d6*=d;
  (*p).d7*=d;
  (*p).d8*=d;
}

void init_exposu3() {
  int k;
  double fctr = 1.0;
  g_exposu3_no_c = 0;
  
  while (fctr>DBL_EPSILON) {
    g_exposu3_no_c++;
    fctr/=(double)(g_exposu3_no_c);
  }
  g_exposu3_no_c += 7;
  g_exposu3_no_c += (g_exposu3_no_c%2);
  
  g_exposu3_c=malloc((g_exposu3_no_c+1)*sizeof(*g_exposu3_c));
  
  g_exposu3_c[0]=1.0;
  for (k=0; k < g_exposu3_no_c; k++)
    g_exposu3_c[k+1]=g_exposu3_c[k]/(double)(k+1);
}

void exposu3(su3* const vr, const su3adj* const p) {
  int n,m,mm;
  su3 ALIGN v,v2,vt;
  su3adj pa;
  double ALIGN d,tc;
  _Complex double t;
  _Complex double ALIGN p0,p1,p2;
  _Complex double ALIGN q0,q1,q2;
  
  _make_su3(v,*p);
  _su3_times_su3(v2,v,v);
  tc = -2.0*(v2.c00 +v2.c11+v2.c22);
  
  pa.d1=(*p).d1;
  pa.d2=(*p).d2;
  pa.d3=(*p).d3;
  pa.d4=(*p).d4;
  pa.d5=(*p).d5;
  pa.d6=(*p).d6;
  pa.d7=(*p).d7;
  pa.d8=(*p).d8;
  
  mm=0;
  while (tc>1.0) {
    mul_su3alg(&pa,0.5);
    tc*=0.5;
    mm+=1;
  }
  
  /* it writes 'p=vec(h_{j,mu})' in matrix form 'v'  */
  _make_su3(v,pa);
  /* calculates v^2 */
  _su3_times_su3(v2,v,v);
  /* t= -tr(X^2)/2*/
  t = -0.5*(v2.c00 +v2.c11+v2.c22);
  /* d= -1i * det(X)*/
  d=-imag_det(&pa);
 /*  printf(" d= %.16f and t=%.16f + 1i %.16f \n",d,creal(t),cimag(t));*/
  
  if(fabs(d)>(1.000001*(1.000002-fabs(t))))
    printf("The norm of X is larger than 1 and N = %d \n", g_exposu3_no_c);
  
  
  p0=g_exposu3_c[g_exposu3_no_c];
  p1=0.0;
  p2=0.0;
  
  for (n=(g_exposu3_no_c-1);n>=0;n--) {
    q0=p0;
    q1=p1;
    q2=p2;
    
    p0=g_exposu3_c[n]-I*d*q2;
    p1=q0-t*q2;
    p2=q1;
  }
   
  /* vr = a0 + a1*v + a2*v2 */
  vt.c00 = p0 + p1 * v.c00 + p2 * v2.c00;
  vt.c01 =      p1 * v.c01 + p2 * v2.c01;
  vt.c02 =      p1 * v.c02 + p2 * v2.c02;
  vt.c10 =      p1 * v.c10 + p2 * v2.c10;
  vt.c11 = p0 + p1 * v.c11 + p2 * v2.c11;
  vt.c12 =      p1 * v.c12 + p2 * v2.c12;
  vt.c20 =      p1 * v.c20 + p2 * v2.c20;
  vt.c21 =      p1 * v.c21 + p2 * v2.c21;
  vt.c22 = p0 + p1 * v.c22 + p2 * v2.c22;
  
  for(m=0;m<mm;m++) {
    _su3_times_su3(v2,vt,vt);
    vt=v2;
  }
  
  vr->c00=vt.c00;
  vr->c01=vt.c01; 
  vr->c02=vt.c02; 
  vr->c10=vt.c10;
  vr->c11=vt.c11;
  vr->c12=vt.c12;
  vr->c20=vt.c20;
  vr->c21=vt.c21;
  vr->c22=vt.c22;
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

  /* compute  row 2 as the conjugate of the cross-product of 3 and 1 */
  vr->c10 = conj(vr->c21 * vr->c02 - vr->c22 * vr->c01);
  vr->c11 = conj(vr->c22 * vr->c00 - vr->c20 * vr->c02);
  vr->c12 = conj(vr->c20 * vr->c01 - vr->c21 * vr->c00);

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

  /* compute  row 2 as the conjugate of the cross-product of 3 and 1 */
  u->c10 = conj(u->c21 * u->c02 - u->c22 * u->c01);
  u->c11 = conj(u->c22 * u->c00 - u->c20 * u->c02);
  u->c12 = conj(u->c20 * u->c01 - u->c21 * u->c00);

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
