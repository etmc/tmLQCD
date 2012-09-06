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
/*******************************************************************************
 *
 *   void diff(spinor * const Q,spinor * const R,spinor * const S)
 *     Makes the difference (*Q) = (*R) - (*S)
 * 
 * Author: Thomas Chiarappa
 *         Thomas.Chiarappa@mib.infn.it
 * 
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "diff_bi.h"

void diff_bi(bispinor * const Q, bispinor * const R, bispinor * const S, const int N){
#ifdef OMP
#pragma omp parallel
  {
#endif

   int ix;
   spinor *q,*r,*s;

/* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
#ifdef OMP
#pragma omp for
#endif
   for (ix = 0; ix < N; ix++) {

     q = (spinor *) &Q[ix].sp_up;
     r = (spinor *) &R[ix].sp_up;
     s = (spinor *) &S[ix].sp_up;
     
     q->s0.c0 = r->s0.c0 - s->s0.c0;
     q->s0.c1 = r->s0.c1 - s->s0.c1;
     q->s0.c2 = r->s0.c2 - s->s0.c2;
     
     q->s1.c0 = r->s1.c0 - s->s1.c0;
     q->s1.c1 = r->s1.c1 - s->s1.c1;
     q->s1.c2 = r->s1.c2 - s->s1.c2;
     
     q->s2.c0 = r->s2.c0 - s->s2.c0;
     q->s2.c1 = r->s2.c1 - s->s2.c1;
     q->s2.c2 = r->s2.c2 - s->s2.c2;
     
     q->s3.c0 = r->s3.c0 - s->s3.c0;
     q->s3.c1 = r->s3.c1 - s->s3.c1;
     q->s3.c2 = r->s3.c2 - s->s3.c2;

     q = (spinor *) &Q[ix].sp_dn;
     r = (spinor *) &R[ix].sp_dn;
     s = (spinor *) &S[ix].sp_dn;
     
     q->s0.c0 = r->s0.c0 - s->s0.c0;
     q->s0.c1 = r->s0.c1 - s->s0.c1;
     q->s0.c2 = r->s0.c2 - s->s0.c2;
     
     q->s1.c0 = r->s1.c0 - s->s1.c0;
     q->s1.c1 = r->s1.c1 - s->s1.c1;
     q->s1.c2 = r->s1.c2 - s->s1.c2;
     
     q->s2.c0 = r->s2.c0 - s->s2.c0;
     q->s2.c1 = r->s2.c1 - s->s2.c1;
     q->s2.c2 = r->s2.c2 - s->s2.c2;
     
     q->s3.c0 = r->s3.c0 - s->s3.c0;
     q->s3.c1 = r->s3.c1 - s->s3.c1;
     q->s3.c2 = r->s3.c2 - s->s3.c2;
   }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
}

