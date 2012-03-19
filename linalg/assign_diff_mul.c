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
#include "su3.h"
#include "assign_diff_mul.h"

/* R=R-c*S */
void assign_diff_mul(spinor * const R, spinor * const S, const _Complex double c, const int N)
{
  spinor *r, *s;

  for (int ix = 0; ix < N; ++ix)
  {
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;

    r->s0.c0 -= c * s->s0.c0;
    r->s0.c1 -= c * s->s0.c1;
    r->s0.c2 -= c * s->s0.c2;

    r->s1.c0 -= c * s->s1.c0;
    r->s1.c1 -= c * s->s1.c1;
    r->s1.c2 -= c * s->s1.c2;

    r->s2.c0 -= c * s->s2.c0;
    r->s2.c1 -= c * s->s2.c1;
    r->s2.c2 -= c * s->s2.c2;

    r->s3.c0 -= c * s->s3.c0;
    r->s3.c1 -= c * s->s3.c1;
    r->s3.c2 -= c * s->s3.c2;
  }
}
