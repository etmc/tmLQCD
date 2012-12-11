/**********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * and modified and extended by Carsten Urbach from 2003-2008
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
 *
 * Hopping_Matrix is the conventional Wilson 
 * hopping matrix
 *
 ****************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#ifdef OMP
#include <omp.h>
#endif
#include <complex.h>
#include "global.h"
#include "su3.h"
#ifdef BGQ
#  include"DirectPut.h"
#endif
#ifdef MPI
#  include "xchange/xchange.h"
#endif
#include "boundary.h"
#include "init/init_dirac_halfspinor.h"
#include "update_backward_gauge.h"
#include "tm_times_Hopping_Matrix.h"

// now comes the definition of tm_times_Hopping_Matrix
// which does (a + g5 i b) * Hopping_Matrix
// where cfactor = a + i b
//

#if (defined _USE_HALFSPINOR && !defined _NO_COMM)
#  include "operator/halfspinor_hopping.h"

#  if ((defined SSE2)||(defined SSE3))
#    include "sse.h"

#  elif (defined BGL && defined XLC)
#    include "bgl.h"

#  elif (defined BGQ && defined XLC)
#    include "bgq.h"
#    include "bgq2.h"
#    include "xlc_prefetch.h"

#  endif

void tm_times_Hopping_Matrix(const int ieo, spinor * const l, spinor * const k, complex double const cfactor) {
  
#  ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#  endif
  
#  ifdef OMP
#  pragma omp parallel
  {
    su3 * restrict u0 ALIGN;
#  endif

#  define _MUL_G5_CMPLX
#  if (defined BGQ && defined XLC)
    complex double ALIGN bla = cfactor;
    vector4double ALIGN cf = vec_ld2(0, (double*) &bla);
#  elif (defined SSE2 || defined SSE3)
    _Complex double ALIGN cf = cfactor;
#  endif
#  include "operator/halfspinor_body.c"
#  undef _MUL_G5_CMPLX    
#  ifdef OMP
  } /* OpenMP closing brace */
#  endif
  return;
}

#elif (!defined _NO_COMM && !defined _USE_HALFSPINOR)
#  include "operator/hopping.h"
#  if ((defined SSE2)||(defined SSE3))
#    include "sse.h"

#  elif (defined BGL && defined XLC)
#    include "bgl.h"

#  elif (defined BGQ && defined XLC)
#    include "bgq.h"
#    include "bgq2.h"
#    include "xlc_prefetch.h"

#  elif defined XLC
#    include"xlc_prefetch.h"

#  endif
void tm_times_Hopping_Matrix(const int ieo, spinor * const l, spinor * const k, double complex const cfactor) {
#  ifdef XLC
#    pragma disjoint(*l, *k)
#  endif
#  ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#  endif

#  if (defined MPI)
  xchange_field(k, ieo);
#  endif
  
#  ifdef OMP
#    pragma omp parallel
  {
#  endif
#  define _MUL_G5_CMPLX
#  if (defined BGQ && defined XLC)
    complex double ALIGN bla = cfactor;
    vector4double ALIGN cf = vec_ld2(0, (double*) &bla);
#  elif (defined SSE2 || defined SSE3)
    _Complex double ALIGN cf = cfactor;
#  endif
#  include "operator/hopping_body_dbl.c"
#  undef _MUL_G5_CMPLX
#  ifdef OMP
  } /* OpenMP closing brace */
#  endif
  return;
}
#endif


