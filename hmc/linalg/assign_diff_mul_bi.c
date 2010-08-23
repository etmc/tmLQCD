/***********************************************************************
 * Copyright (C) 2006 Thomas Chiarappa
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
/* $Id$ */

/************************************************************************
 * 
 *      Adpated routine evaluating the S=S-c*Q where S,Q are bispinors
 * 
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include "su3.h"
#include "assign_diff_mul_bi.h"


/* S=S-c*Q */
void assign_diff_mul_bi(bispinor * const S, bispinor * const R, const complex c, const int N){


  int ix;
  double a,b;
  spinor *r, *s;

  a = - c.re;
  b = - c.im;

  for (ix=0;ix<N;ix++){

    s = (spinor *) &S[ix].sp_up;
    r = (spinor *) &R[ix].sp_up;

    (*s).s0.c0.re += a*(*r).s0.c0.re - b*(*r).s0.c0.im;
    (*s).s0.c0.im += a*(*r).s0.c0.im + b*(*r).s0.c0.re;
    (*s).s0.c1.re += a*(*r).s0.c1.re - b*(*r).s0.c1.im;
    (*s).s0.c1.im += a*(*r).s0.c1.im + b*(*r).s0.c1.re;
    (*s).s0.c2.re += a*(*r).s0.c2.re - b*(*r).s0.c2.im;
    (*s).s0.c2.im += a*(*r).s0.c2.im + b*(*r).s0.c2.re;

    (*s).s1.c0.re += a*(*r).s1.c0.re - b*(*r).s1.c0.im;
    (*s).s1.c0.im += a*(*r).s1.c0.im + b*(*r).s1.c0.re;
    (*s).s1.c1.re += a*(*r).s1.c1.re - b*(*r).s1.c1.im;
    (*s).s1.c1.im += a*(*r).s1.c1.im + b*(*r).s1.c1.re;
    (*s).s1.c2.re += a*(*r).s1.c2.re - b*(*r).s1.c2.im;
    (*s).s1.c2.im += a*(*r).s1.c2.im + b*(*r).s1.c2.re;

    (*s).s2.c0.re += a*(*r).s2.c0.re - b*(*r).s2.c0.im;
    (*s).s2.c0.im += a*(*r).s2.c0.im + b*(*r).s2.c0.re;
    (*s).s2.c1.re += a*(*r).s2.c1.re - b*(*r).s2.c1.im;
    (*s).s2.c1.im += a*(*r).s2.c1.im + b*(*r).s2.c1.re;
    (*s).s2.c2.re += a*(*r).s2.c2.re - b*(*r).s2.c2.im;
    (*s).s2.c2.im += a*(*r).s2.c2.im + b*(*r).s2.c2.re;
     
    (*s).s3.c0.re += a*(*r).s3.c0.re - b*(*r).s3.c0.im;
    (*s).s3.c0.im += a*(*r).s3.c0.im + b*(*r).s3.c0.re;
    (*s).s3.c1.re += a*(*r).s3.c1.re - b*(*r).s3.c1.im;
    (*s).s3.c1.im += a*(*r).s3.c1.im + b*(*r).s3.c1.re;
    (*s).s3.c2.re += a*(*r).s3.c2.re - b*(*r).s3.c2.im;
    (*s).s3.c2.im += a*(*r).s3.c2.im + b*(*r).s3.c2.re;


    s = (spinor *) &S[ix].sp_dn;
    r = (spinor *) &R[ix].sp_dn;

    (*s).s0.c0.re += a*(*r).s0.c0.re - b*(*r).s0.c0.im;
    (*s).s0.c0.im += a*(*r).s0.c0.im + b*(*r).s0.c0.re;
    (*s).s0.c1.re += a*(*r).s0.c1.re - b*(*r).s0.c1.im;
    (*s).s0.c1.im += a*(*r).s0.c1.im + b*(*r).s0.c1.re;
    (*s).s0.c2.re += a*(*r).s0.c2.re - b*(*r).s0.c2.im;
    (*s).s0.c2.im += a*(*r).s0.c2.im + b*(*r).s0.c2.re;

    (*s).s1.c0.re += a*(*r).s1.c0.re - b*(*r).s1.c0.im;
    (*s).s1.c0.im += a*(*r).s1.c0.im + b*(*r).s1.c0.re;
    (*s).s1.c1.re += a*(*r).s1.c1.re - b*(*r).s1.c1.im;
    (*s).s1.c1.im += a*(*r).s1.c1.im + b*(*r).s1.c1.re;
    (*s).s1.c2.re += a*(*r).s1.c2.re - b*(*r).s1.c2.im;
    (*s).s1.c2.im += a*(*r).s1.c2.im + b*(*r).s1.c2.re;

    (*s).s2.c0.re += a*(*r).s2.c0.re - b*(*r).s2.c0.im;
    (*s).s2.c0.im += a*(*r).s2.c0.im + b*(*r).s2.c0.re;
    (*s).s2.c1.re += a*(*r).s2.c1.re - b*(*r).s2.c1.im;
    (*s).s2.c1.im += a*(*r).s2.c1.im + b*(*r).s2.c1.re;
    (*s).s2.c2.re += a*(*r).s2.c2.re - b*(*r).s2.c2.im;
    (*s).s2.c2.im += a*(*r).s2.c2.im + b*(*r).s2.c2.re;
     
    (*s).s3.c0.re += a*(*r).s3.c0.re - b*(*r).s3.c0.im;
    (*s).s3.c0.im += a*(*r).s3.c0.im + b*(*r).s3.c0.re;
    (*s).s3.c1.re += a*(*r).s3.c1.re - b*(*r).s3.c1.im;
    (*s).s3.c1.im += a*(*r).s3.c1.im + b*(*r).s3.c1.re;
    (*s).s3.c2.re += a*(*r).s3.c2.re - b*(*r).s3.c2.im;
    (*s).s3.c2.im += a*(*r).s3.c2.im + b*(*r).s3.c2.re;

  }
}


