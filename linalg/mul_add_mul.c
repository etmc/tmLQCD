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
#include "su3.h"
#include "mul_add_mul.h"


/* Makes (*R)=c1*(*S)+c2*(*U) , c1 and c2 are complex constants */
void mul_add_mul(spinor * const R,spinor * const S,spinor * const U,const _Complex double c1,const _Complex double c2, const int N)
{
  spinor *r, *s, *u;
  
  for (int ix = 0; ix < N; ++ix)
  {
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    u=(spinor *) U + ix;
    
    r->s0.c0 = c1 * s->s0.c0 + c2 * u->s0.c0;
    r->s0.c1 = c1 * s->s0.c1 + c2 * u->s0.c1;
    r->s0.c2 = c1 * s->s0.c2 + c2 * u->s0.c2;

    r->s1.c0 = c1 * s->s1.c0 + c2 * u->s1.c0;
    r->s1.c1 = c1 * s->s1.c1 + c2 * u->s1.c1;
    r->s1.c2 = c1 * s->s1.c2 + c2 * u->s1.c2;
    
    r->s2.c0 = c1 * s->s2.c0 + c2 * u->s2.c0;
    r->s2.c1 = c1 * s->s2.c1 + c2 * u->s2.c1;
    r->s2.c2 = c1 * s->s2.c2 + c2 * u->s2.c2;
    
    r->s3.c0 = c1 * s->s3.c0 + c2 * u->s3.c0;
    r->s3.c1 = c1 * s->s3.c1 + c2 * u->s3.c1;
    r->s3.c2 = c1 * s->s3.c2 + c2 * u->s3.c2;
  }
}
