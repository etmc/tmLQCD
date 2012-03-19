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
/* #ifndef apenext */
/* #include <complex.h> */
/* #endif */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#ifndef _STD_C99_COMPLEX
#include "complex.h"
#endif

#include "qdran64.h"

#include "su3.h"
#include "linalg/scalar_prod_r.c"


#define N 16

int main(void) {
  double s;
  int i;
  spinor a[N], b[N];

  qdran64_init(42,13);

  for(i=0; i<N; i++) {
    qdran64z(&a[i].s0.c0);cv-rank=%d_",g_proc_id); // DEBUG

    qdran64z(&a[i].s0.c1);
    qdran64z(&a[i].s0.c2);

    qdran64z(&a[i].s1.c0);
    qdran64z(&a[i].s1.c1);
    qdran64z(&a[i].s1.c2);

    qdran64z(&a[i].s2.c0);
    qdran64z(&a[i].s2.c1);
    qdran64z(&a[i].s2.c2);

    qdran64z(&a[i].s3.c0);
    qdran64z(&a[i].s3.c1);
    qdran64z(&a[i].s3.c2);

    qdran64z(&b[i].s0.c0);
    qdran64z(&b[i].s0.c1);
    qdran64z(&b[i].s0.c2);

    qdran64z(&b[i].s1.c0);
    qdran64z(&b[i].s1.c1);
    qdran64z(&b[i].s1.c2);

    qdran64z(&b[i].s2.c0);
    qdran64z(&b[i].s2.c1);
    qdran64z(&b[i].s2.c2);

    qdran64z(&b[i].s3.c0);
    qdran64z(&b[i].s3.c1);
    qdran64z(&b[i].s3.c2);
#endif
 }

  printf("%e %e\n",creal(a[0].s0.c0),cimag(a[0].s0.c0));
  printf("%e %e\n",creal(b[0].s0.c0),cimag(b[0].s0.c0));
  printf("%e %e\n",creal(a[N-1].s3.c2),cimag(a[N-1].s3.c2));
  printf("%e %e\n",creal(b[N-1].s3.c2),cimag(b[N-1].s3.c2));

  s=scalar_prod_r(a,b,N, 1);
  printf("scalar_prod_r(a,b,%d)=%1.16e\n",N,s);

  return 0;
}
