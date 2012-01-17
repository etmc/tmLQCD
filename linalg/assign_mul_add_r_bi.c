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
/* $Id: assign_mul_add_r_bi.c 1150 2009-02-16 16:52:09Z urbach $ */


/************************************************************************
 *
 *      Adpated routine evaluating the P=c*P+Q where P,Q are bispinors
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
#include "sse.h"
#include "assign_mul_add_r_bi.h"


#if defined SSE2
/* k input , l output*/
void assign_mul_add_r_bi(bispinor * const S, const double c, bispinor * const R, const int N) {

  int ix;
  su3_vector *s,*r;
  /*
  su3_vector *t,*u;
  */

  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c));
  
  
  s=(su3_vector *) &S[0].sp_up.s0;
  r=(su3_vector *) &R[0].sp_up.s0;

/*  for (ix=0;ix<4*N;ix++) { */
  for (ix=0;ix<2*4*N;ix++) {
    _sse_load(*s);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm0 \n\t"
			  "mulpd %%xmm7, %%xmm1 \n\t"
			  "mulpd %%xmm7, %%xmm2"
			  :
			  :);
    _sse_load_up(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*s);
    s++; r++;
  }

}

#else
/* k input , l output*/
void assign_mul_add_r_bi(bispinor * const R, const double c, bispinor * const S, const int N)
{
  spinor *r,*s;
  
  /* Change due to even-odd preconditioning : VOLUME   to VOLUME/2 */   
  for (int ix = 0; ix < N; ++ix)
  {
    r = (spinor *) &R[ix].sp_up;
    s = (spinor *) &S[ix].sp_up;
    
    r->s0.c0 = c * r->s0.c0 + s->s0.c0;
    r->s0.c1 = c * r->s0.c1 + s->s0.c1;
    r->s0.c2 = c * r->s0.c2 + s->s0.c2;    

    r->s1.c0 = c * r->s1.c0 + s->s1.c0;
    r->s1.c1 = c * r->s1.c1 + s->s1.c1;
    r->s1.c2 = c * r->s1.c2 + s->s1.c2;    

    r->s2.c0 = c * r->s2.c0 + s->s2.c0;
    r->s2.c1 = c * r->s2.c1 + s->s2.c1;
    r->s2.c2 = c * r->s2.c2 + s->s2.c2;    

    r->s3.c0 = c * r->s3.c0 + s->s3.c0;
    r->s3.c1 = c * r->s3.c1 + s->s3.c1;
    r->s3.c2 = c * r->s3.c2 + s->s3.c2;    

    r = (spinor *) &R[ix].sp_dn;
    s = (spinor *) &S[ix].sp_dn;
    
    r->s0.c0 = c * r->s0.c0 + s->s0.c0;
    r->s0.c1 = c * r->s0.c1 + s->s0.c1;
    r->s0.c2 = c * r->s0.c2 + s->s0.c2;    

    r->s1.c0 = c * r->s1.c0 + s->s1.c0;
    r->s1.c1 = c * r->s1.c1 + s->s1.c1;
    r->s1.c2 = c * r->s1.c2 + s->s1.c2;    

    r->s2.c0 = c * r->s2.c0 + s->s2.c0;
    r->s2.c1 = c * r->s2.c1 + s->s2.c1;
    r->s2.c2 = c * r->s2.c2 + s->s2.c2;    

    r->s3.c0 = c * r->s3.c0 + s->s3.c0;
    r->s3.c1 = c * r->s3.c1 + s->s3.c1;
    r->s3.c2 = c * r->s3.c2 + s->s3.c2;    
  }
}

#endif



