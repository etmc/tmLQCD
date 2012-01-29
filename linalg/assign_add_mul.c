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
 * File assign_add_mul.c 
 *
 *   void assign_add_mul(spinor * const P, spinor * const Q, const complex c)
 *     (*P) = (*P) + c(*Q)        c is a complex constant
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_add_mul.h"


void assign_add_mul(spinor * const P, spinor * const Q, const complex c, const int N){
  int ix;
  static double a,b;
  spinor *r,*s;

  a=c.re;
  b=c.im;
   
  for (ix=0;ix<N;ix++){
    r=(spinor *) P + ix;
    s=(spinor *) Q + ix;

    (*r).s0.c0.re+=a*(*s).s0.c0.re - b*(*s).s0.c0.im;
    (*r).s0.c0.im+=a*(*s).s0.c0.im + b*(*s).s0.c0.re;
    (*r).s0.c1.re+=a*(*s).s0.c1.re - b*(*s).s0.c1.im;
    (*r).s0.c1.im+=a*(*s).s0.c1.im + b*(*s).s0.c1.re;
    (*r).s0.c2.re+=a*(*s).s0.c2.re - b*(*s).s0.c2.im;
    (*r).s0.c2.im+=a*(*s).s0.c2.im + b*(*s).s0.c2.re;

    (*r).s1.c0.re+=a*(*s).s1.c0.re - b*(*s).s1.c0.im;
    (*r).s1.c0.im+=a*(*s).s1.c0.im + b*(*s).s1.c0.re;
    (*r).s1.c1.re+=a*(*s).s1.c1.re - b*(*s).s1.c1.im;
    (*r).s1.c1.im+=a*(*s).s1.c1.im + b*(*s).s1.c1.re;
    (*r).s1.c2.re+=a*(*s).s1.c2.re - b*(*s).s1.c2.im;
    (*r).s1.c2.im+=a*(*s).s1.c2.im + b*(*s).s1.c2.re;

    (*r).s2.c0.re+=a*(*s).s2.c0.re - b*(*s).s2.c0.im;
    (*r).s2.c0.im+=a*(*s).s2.c0.im + b*(*s).s2.c0.re;
    (*r).s2.c1.re+=a*(*s).s2.c1.re - b*(*s).s2.c1.im;
    (*r).s2.c1.im+=a*(*s).s2.c1.im + b*(*s).s2.c1.re;
    (*r).s2.c2.re+=a*(*s).s2.c2.re - b*(*s).s2.c2.im;
    (*r).s2.c2.im+=a*(*s).s2.c2.im + b*(*s).s2.c2.re;

    (*r).s3.c0.re+=a*(*s).s3.c0.re - b*(*s).s3.c0.im;
    (*r).s3.c0.im+=a*(*s).s3.c0.im + b*(*s).s3.c0.re;
    (*r).s3.c1.re+=a*(*s).s3.c1.re - b*(*s).s3.c1.im;
    (*r).s3.c1.im+=a*(*s).s3.c1.im + b*(*s).s3.c1.re;
    (*r).s3.c2.re+=a*(*s).s3.c2.re - b*(*s).s3.c2.im;
    (*r).s3.c2.im+=a*(*s).s3.c2.im + b*(*s).s3.c2.re;

  }
}


