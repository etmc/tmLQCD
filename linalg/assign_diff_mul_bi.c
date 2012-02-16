/***********************************************************************
 * Copyright (C) 2006 Thomas Chiarappa
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is e softw: you candistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the e Softw Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for m details.
 * 
 * You should haveceived a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

/************************************************************************
 * 
 *      Adapted routine evaluating the S=S-c*Q wh S,Q  bispinors
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
void assign_diff_mul_bi(bispinor * const S, bispinor * const R, const _Complex double c, const int N)
{
  spinor *r, *s;

  for (int ix = 0; ix < N; ++ix)
  {
    s = (spinor *) &S[ix].sp_up;
    r = (spinor *) &R[ix].sp_up;

    s->s0.c0 -= c * r->s0.c0;
    s->s0.c1 -= c * r->s0.c1;
    s->s0.c2 -= c * r->s0.c2;

    s->s1.c0 -= c * r->s1.c0;
    s->s1.c1 -= c * r->s1.c1;
    s->s1.c2 -= c * r->s1.c2;

    s->s2.c0 -= c * r->s2.c0;
    s->s2.c1 -= c * r->s2.c1;
    s->s2.c2 -= c * r->s2.c2;
     
    s->s3.c0 -= c * r->s3.c0;
    s->s3.c1 -= c * r->s3.c1;
    s->s3.c2 -= c * r->s3.c2;


    s = (spinor *) &S[ix].sp_dn;
    r = (spinor *) &R[ix].sp_dn;

    s->s0.c0 -= c * r->s0.c0;
    s->s0.c1 -= c * r->s0.c1;
    s->s0.c2 -= c * r->s0.c2;

    s->s1.c0 -= c * r->s1.c0;
    s->s1.c1 -= c * r->s1.c1;
    s->s1.c2 -= c * r->s1.c2;

    s->s2.c0 -= c * r->s2.c0;
    s->s2.c1 -= c * r->s2.c1;
    s->s2.c2 -= c * r->s2.c2;
     
    s->s3.c0 -= c * r->s3.c0;
    s->s3.c1 -= c * r->s3.c1;
    s->s3.c2 -= c * r->s3.c2;
  }
}


