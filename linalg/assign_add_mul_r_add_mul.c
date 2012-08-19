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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "su3.h"
#include "su3adj.h"
#include "sse.h"
#include "assign_add_mul_r_add_mul.h"

#if ( defined SSE2 || defined SSE3 )
void assign_add_mul_r_add_mul(spinor * const R, spinor * const S, spinor * const U,
			      const double c1,const double c2, const int N) {

  int ix;
  su3_vector *s,*r,*t;
  r=&R[0].s0;
  s=&S[0].s0;
  t=&U[0].s0;
  __asm__ __volatile__ ("movsd %0, %%xmm6 \n\t"
			"unpcklpd %%xmm6, %%xmm6"
			:
			:
			"m" (c1));
  __asm__ __volatile__ ("movsd %0, %%xmm7 \n\t"
			"unpcklpd %%xmm7, %%xmm7"
			:
			:
			"m" (c2));
  
  for (ix = 0; ix < 4*N; ix++) {
    _sse_load_up(*s);
    __asm__ __volatile__ ("mulpd %%xmm6, %%xmm3 \n\t"
			  "mulpd %%xmm6, %%xmm4 \n\t"
			  "mulpd %%xmm6, %%xmm5"
			  :
			  :);
    _sse_load(*r);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_load_up(*t);
    __asm__ __volatile__ ("mulpd %%xmm7, %%xmm3 \n\t"
			  "mulpd %%xmm7, %%xmm4 \n\t"
			  "mulpd %%xmm7, %%xmm5"
			  :
			  :);
    __asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t"
			  "addpd %%xmm4, %%xmm1 \n\t"
			  "addpd %%xmm5, %%xmm2"
			  :
			  :);
    _sse_store(*r);
    r++; s++; t++;
  }
}
#else
/* j, k input, l output */
void assign_add_mul_r_add_mul(spinor * const R, spinor * const S, spinor * const U,
			      const double c1,const double c2, const int N) {
#ifdef OMP
#pragma omp parallel
  {
#endif

   spinor *r,*s,*t;

#ifdef OMP
#pragma omp for
#endif
   for (int ix = 0; ix < N; ++ix)
   {
     r=R+ix;      
     s=S+ix;
     t=U+ix;
     
     r->s0.c0 += c1 * s->s0.c0 + c2 * t->s0.c0;
     r->s0.c1 += c1 * s->s0.c1 + c2 * t->s0.c1;
     r->s0.c2 += c1 * s->s0.c2 + c2 * t->s0.c2;

     r->s1.c0 += c1 * s->s1.c0 + c2 * t->s1.c0;
     r->s1.c1 += c1 * s->s1.c1 + c2 * t->s1.c1;
     r->s1.c2 += c1 * s->s1.c2 + c2 * t->s1.c2;

     r->s2.c0 += c1 * s->s2.c0 + c2 * t->s2.c0;
     r->s2.c1 += c1 * s->s2.c1 + c2 * t->s2.c1;
     r->s2.c2 += c1 * s->s2.c2 + c2 * t->s2.c2;
     
     r->s3.c0 += c1 * s->s3.c0 + c2 * t->s3.c0;
     r->s3.c1 += c1 * s->s3.c1 + c2 * t->s3.c1;
     r->s3.c2 += c1 * s->s3.c2 + c2 * t->s3.c2;
   }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

}
#endif
