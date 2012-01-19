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
/*******************************************************************************
 *
 * File assign_add_mul_add_mul.c
 *
 *   void assign_add_mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2)
 *     (*R) = (*R) + c1*(*S) + c2*(*U) with c1 and c2 complex variables
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_add_mul_add_mul.h"


/* S,U input, R inoutput, c1,c2 input */
void assign_add_mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const complex c1,const complex c2, const int N){
  int ix;
  spinor *r,*s,*u;
  double a1,b1,a2,b2;

  a1 = c1.re;
  b1 = c1.im;

  a2 = c2.re;
  b2 = c2.im;

  for (ix=0;ix<N;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    u=(spinor *) U + ix;
    
    (*r).s0.c0.re+=a1*(*s).s0.c0.re - b1*(*s).s0.c0.im + a2*(*u).s0.c0.re - b2*(*u).s0.c0.im;
    (*r).s0.c1.re+=a1*(*s).s0.c1.re - b1*(*s).s0.c1.im + a2*(*u).s0.c1.re - b2*(*u).s0.c1.im;
    (*r).s0.c2.re+=a1*(*s).s0.c2.re - b1*(*s).s0.c2.im + a2*(*u).s0.c2.re - b2*(*u).s0.c2.im;

    (*r).s1.c0.re+=a1*(*s).s1.c0.re - b1*(*s).s1.c0.im + a2*(*u).s1.c0.re - b2*(*u).s1.c0.im;
    (*r).s1.c1.re+=a1*(*s).s1.c1.re - b1*(*s).s1.c1.im + a2*(*u).s1.c1.re - b2*(*u).s1.c1.im;
    (*r).s1.c2.re+=a1*(*s).s1.c2.re - b1*(*s).s1.c2.im + a2*(*u).s1.c2.re - b2*(*u).s1.c2.im;

    (*r).s2.c0.re+=a1*(*s).s2.c0.re - b1*(*s).s2.c0.im + a2*(*u).s2.c0.re - b2*(*u).s2.c0.im;
    (*r).s2.c1.re+=a1*(*s).s2.c1.re - b1*(*s).s2.c1.im + a2*(*u).s2.c1.re - b2*(*u).s2.c1.im;
    (*r).s2.c2.re+=a1*(*s).s2.c2.re - b1*(*s).s2.c2.im + a2*(*u).s2.c2.re - b2*(*u).s2.c2.im;

    (*r).s3.c0.re+=a1*(*s).s3.c0.re - b1*(*s).s3.c0.im + a2*(*u).s3.c0.re - b2*(*u).s3.c0.im;
    (*r).s3.c1.re+=a1*(*s).s3.c1.re - b1*(*s).s3.c1.im + a2*(*u).s3.c1.re - b2*(*u).s3.c1.im;
    (*r).s3.c2.re+=a1*(*s).s3.c2.re - b1*(*s).s3.c2.im + a2*(*u).s3.c2.re - b2*(*u).s3.c2.im;



    (*r).s0.c0.im+=a1*(*s).s0.c0.im + b1*(*s).s0.c0.re + a2*(*u).s0.c0.im + b2*(*u).s0.c0.re;
    (*r).s0.c1.im+=a1*(*s).s0.c1.im + b1*(*s).s0.c1.re + a2*(*u).s0.c1.im + b2*(*u).s0.c1.re;
    (*r).s0.c2.im+=a1*(*s).s0.c2.im + b1*(*s).s0.c2.re + a2*(*u).s0.c2.im + b2*(*u).s0.c2.re;

    (*r).s1.c0.im+=a1*(*s).s1.c0.im + b1*(*s).s1.c0.re + a2*(*u).s1.c0.im + b2*(*u).s1.c0.re;
    (*r).s1.c1.im+=a1*(*s).s1.c1.im + b1*(*s).s1.c1.re + a2*(*u).s1.c1.im + b2*(*u).s1.c1.re;
    (*r).s1.c2.im+=a1*(*s).s1.c2.im + b1*(*s).s1.c2.re + a2*(*u).s1.c2.im + b2*(*u).s1.c2.re;

    (*r).s2.c0.im+=a1*(*s).s2.c0.im + b1*(*s).s2.c0.re + a2*(*u).s2.c0.im + b2*(*u).s2.c0.re;
    (*r).s2.c1.im+=a1*(*s).s2.c1.im + b1*(*s).s2.c1.re + a2*(*u).s2.c1.im + b2*(*u).s2.c1.re;
    (*r).s2.c2.im+=a1*(*s).s2.c2.im + b1*(*s).s2.c2.re + a2*(*u).s2.c2.im + b2*(*u).s2.c2.re;

    (*r).s3.c0.im+=a1*(*s).s3.c0.im + b1*(*s).s3.c0.re + a2*(*u).s3.c0.im + b2*(*u).s3.c0.re;
    (*r).s3.c1.im+=a1*(*s).s3.c1.im + b1*(*s).s3.c1.re + a2*(*u).s3.c1.im + b2*(*u).s3.c1.re;
    (*r).s3.c2.im+=a1*(*s).s3.c2.im + b1*(*s).s3.c2.re + a2*(*u).s3.c2.im + b2*(*u).s3.c2.re;

  }
}
