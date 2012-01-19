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
 * File assign.c
 *
 *   void assign(spinor * const R, spinor * const S)
 *     Assign (*R) = (*S)
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#ifdef _STD_C99_COMPLEX_CHECKED
# include <complex.h>
#endif
#ifdef apenext
# include <topology.h>
# include <queue.h>
#endif
#include <math.h>
#include "su3.h"
#include "assign.h"



/* S input, R output */

#if ((!defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void assign(spinor * const R, spinor * const S, const int N){
  int ix;
  spinor *r,*s;
  
  for (ix = 0; ix < N; ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s0.c0.re = (*s).s0.c0.re;
    (*r).s0.c0.im = (*s).s0.c0.im;
    (*r).s0.c1.re = (*s).s0.c1.re;
    (*r).s0.c1.im = (*s).s0.c1.im;
    (*r).s0.c2.re = (*s).s0.c2.re;
    (*r).s0.c2.im = (*s).s0.c2.im;
    
    (*r).s1.c0.re = (*s).s1.c0.re;
    (*r).s1.c0.im = (*s).s1.c0.im;
    (*r).s1.c1.re = (*s).s1.c1.re;
    (*r).s1.c1.im = (*s).s1.c1.im;
    (*r).s1.c2.re = (*s).s1.c2.re;
    (*r).s1.c2.im = (*s).s1.c2.im;         
    
    (*r).s2.c0.re = (*s).s2.c0.re;
    (*r).s2.c0.im = (*s).s2.c0.im;
    (*r).s2.c1.re = (*s).s2.c1.re;
    (*r).s2.c1.im = (*s).s2.c1.im;
    (*r).s2.c2.re = (*s).s2.c2.re;
    (*r).s2.c2.im = (*s).s2.c2.im;         
    
    (*r).s3.c0.re = (*s).s3.c0.re;
    (*r).s3.c0.im = (*s).s3.c0.im;
    (*r).s3.c1.re = (*s).s3.c1.re;
    (*r).s3.c1.im = (*s).s3.c1.im;
    (*r).s3.c2.re = (*s).s3.c2.re;
    (*r).s3.c2.im = (*s).s3.c2.im;
  }
}
#endif

#if ((defined _STD_C99_COMPLEX_CHECKED) && (!defined apenext))

void assign(spinor * const R, spinor * const S, const int N){
  register int ix=0;
  register spinor *rPointer,*sPointer;

  rPointer = R;
  sPointer = S;

  do {
    register spinor s, r;
    ix+=1;

    s = *(sPointer);

    r.s0.c0 = s.s0.c0;
    r.s0.c1 = s.s0.c1;
    r.s0.c2 = s.s0.c2;

    r.s1.c0 = s.s1.c0;
    r.s1.c1 = s.s1.c1;
    r.s1.c2 = s.s1.c2;

    r.s2.c0 = s.s2.c0;
    r.s2.c1 = s.s2.c1;
    r.s2.c2 = s.s2.c2;
    
    r.s3.c0 = s.s3.c0;
    r.s3.c1 = s.s3.c1;
    r.s3.c2 = s.s3.c2;

    *(rPointer) = r;

    rPointer+=1;
    sPointer+=1;

 } while (ix<N);
}
#endif


#ifdef apenext

#define NOWHERE_COND(condition) ((condition) ? 0x0 : NOWHERE )

void assign(spinor * const R, spinor * const S, const int N){
  register int ix=N;
  register spinor *rPointer,*sPointer;

  rPointer = R;
  sPointer = S;

  prefetch (*(sPointer));
  {
#pragma cache

    do {
      register spinor s, r;
      ix-=1;

      sPointer+=1;

      fetch(s);

      prefetch (*(sPointer+NOWHERE_COND(ix)));
      
      r.s0.c0 = s.s0.c0;
      r.s0.c1 = s.s0.c1;
      r.s0.c2 = s.s0.c2;

      r.s1.c0 = s.s1.c0;
      r.s1.c1 = s.s1.c1;
      r.s1.c2 = s.s1.c2;

      r.s2.c0 = s.s2.c0;
      r.s2.c1 = s.s2.c1;
      r.s2.c2 = s.s2.c2;
    
      r.s3.c0 = s.s3.c0;
      r.s3.c1 = s.s3.c1;
      r.s3.c2 = s.s3.c2;

      *(rPointer) = r;

      rPointer+=1;

    } while (ix>0);
  }
}
#endif

#ifdef WITHLAPH
void assign_su3vect(su3_vector * const R, su3_vector * const S, const int N)
{
int ix;
su3_vector *r,*s;

	for (ix = 0; ix < N; ix++)
	{
    r=(su3_vector *) R + ix;
    s=(su3_vector *) S + ix;
    
    (*r).c0.re = (*s).c0.re;
    (*r).c0.im = (*s).c0.im;
    (*r).c1.re = (*s).c1.re;
    (*r).c1.im = (*s).c1.im;
    (*r).c2.re = (*s).c2.re;
    (*r).c2.im = (*s).c2.im;
	}
}
#endif
