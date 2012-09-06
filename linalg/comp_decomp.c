/***********************************************************************
 *
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
 *
 * File comp_decomp.c
 *
 *   void compact(spinor * const R, spinor * const S, spinor * const T); 
 *     
 *    Builds the Bi-spinor R out of the spinors S and T 
 *        S in first half (top)    T in second half (bottom)
 *
 *
 *   void decompact(spinor * const S, spinor * const T, spinor * const R); 
 *     
 *            Splits the Bi-spinor R in the spinors S and T 
 *        S in first half (top)    T in second half (bottom)
 *
 *****************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "comp_decomp.h"


/* S and P inputs, R output */
void compact(bispinor * const R, spinor * const S, spinor * const P)
{ 
#ifdef OMP
#pragma omp parallel
  {
#endif
  spinor *r,*s;
  spinor *u,*t;
  
#ifdef OMP
#pragma omp for
#endif
  for (int ix = 0; ix < VOLUME/2; ix++){
    r=(spinor *) &R[ix].sp_up;
    s=(spinor *) S + ix;
    
    /*    (*r) = (*s);    */

    r->s0.c0 = s->s0.c0;
    r->s0.c1 = s->s0.c1;
    r->s0.c2 = s->s0.c2;
    
    r->s1.c0 = s->s1.c0;
    r->s1.c1 = s->s1.c1;
    r->s1.c2 = s->s1.c2;
    
    r->s2.c0 = s->s2.c0;
    r->s2.c1 = s->s2.c1;
    r->s2.c2 = s->s2.c2;
    
    r->s3.c0 = s->s3.c0;
    r->s3.c1 = s->s3.c1;
    r->s3.c2 = s->s3.c2;

    u=(spinor *) &R[ix].sp_dn;
    t=(spinor *) P + ix;
    
    u->s0.c0 = t->s0.c0;
    u->s0.c1 = t->s0.c1;
    u->s0.c2 = t->s0.c2;
    
    u->s1.c0 = t->s1.c0;
    u->s1.c1 = t->s1.c1;
    u->s1.c2 = t->s1.c2;
    
    u->s2.c0 = t->s2.c0;
    u->s2.c1 = t->s2.c1;
    u->s2.c2 = t->s2.c2;
    
    u->s3.c0 = t->s3.c0;
    u->s3.c1 = t->s3.c1;
    u->s3.c2 = t->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
  
  /* 
      The following IS NOT enough, since it copies the values 
      starting from the adress given by the pointer &R->sp_up, but
      following the colour - spin - FLAVOUR - volume structure.

      In other words: staring with the FIRST site on the lattice ([0]) 
      the routine copies the first 3 (colour) * 4 (spin) component of 
      the spinor  S  onto the corresponding adresses of the spinor 
      R->sp_up . Then it continues by copying the component 
      S[1].s0.c0  onto the adress  R[0].sp_dn (.s0.c0),  
      !!!  AND NOT JUMPING TO  R[1].sp_up (.s0.c0)  !!!
      because of the structure and mem. allocation of the bispinor
  */

  /*
  assign(&R->sp_up, &S[0], VOLUME/2);
  */
  /*  
  assign(&R->sp_dn, &P[0], VOLUME/2);
  */

}


/* R input , S and P outputs */
void decompact(spinor * const S, spinor * const P, bispinor * const R){
#ifdef OMP
#pragma omp parallel
  {
#endif

  spinor *r,*s;
  spinor *u,*t;
  
#ifdef OMP
#pragma omp for
#endif
  for (int ix = 0; ix < VOLUME/2; ix++)
  {
    s=(spinor *) &R[ix].sp_up;
    r=(spinor *) S + ix;
    
    r->s0.c0 = s->s0.c0;
    r->s0.c1 = s->s0.c1;
    r->s0.c2 = s->s0.c2;
    
    r->s1.c0 = s->s1.c0;
    r->s1.c1 = s->s1.c1;
    r->s1.c2 = s->s1.c2;
    
    r->s2.c0 = s->s2.c0;
    r->s2.c1 = s->s2.c1;
    r->s2.c2 = s->s2.c2;
    
    r->s3.c0 = s->s3.c0;
    r->s3.c1 = s->s3.c1;
    r->s3.c2 = s->s3.c2;


    t=(spinor *) &R[ix].sp_dn;
    u=(spinor *) P + ix;
    
    u->s0.c0 = t->s0.c0;
    u->s0.c1 = t->s0.c1;
    u->s0.c2 = t->s0.c2;
    
    u->s1.c0 = t->s1.c0;
    u->s1.c1 = t->s1.c1;
    u->s1.c2 = t->s1.c2;
    
    u->s2.c0 = t->s2.c0;
    u->s2.c1 = t->s2.c1;
    u->s2.c2 = t->s2.c2;
    
    u->s3.c0 = t->s3.c0;
    u->s3.c1 = t->s3.c1;
    u->s3.c2 = t->s3.c2;
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif

  /* !!! The following should be enough !!! */

 /* 
      The following IS NOT enough,  See explanation above 
 */
  /*
  assign(&S[0], &R->sp_up, VOLUME/2);

  assign(&T[0], &R->sp_dn, VOLUME/2);
  */

}
