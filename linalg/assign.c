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
#include <complex.h>
#include <math.h>
#include "su3.h"
#include "assign.h"

/* S input, R output */
void assign(spinor * const R, spinor * const S, const int N)
{
  memmove(R, S, N * sizeof(spinor));
}

#ifdef WITHLAPH
void assign_su3vect(su3_vector * const R, su3_vector * const S, const int N)
{
  su3_vector *r,*s;

  for (int ix = 0; ix < N; ++ix) 
  {
    r=P+ix;      
    s=Q+ix;
    
    r->c0 = s->c0;
    r->c1 = s->c1;
    r->c2 = s->c2;
  }
}
}
#endif
