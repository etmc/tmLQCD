/**********************************************************************
 *
 *
 * Copyright (C) 2001 Martin Luescher
 *               2002 Martin Hasenbusch
 *               2003, 2004, 2005, 2006, 2007, 2008 Carsten Urbach
 *
 * BG and halfspinor versions (C) 2007, 2008 Carsten Urbach
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
 * \kappa\sum_{\pm\mu}(r+\gamma_\mu)U_{x,\mu}
 *
 * for ieo = 0 this is M_{eo}, for ieo = 1
 * it is M_{oe}
 *
 * l is the output, k the input field
 *
 *  Structure of top level precompiler directives 
 *
 * - defining _USE_HALFSPINOR implies that we also use
 *   a "gauge copy"
 *
 * - such that we are checking for the _USE_GAUGECOPY feature seperatly in the 
 *   ELSE branch of the "if defined _USE_HALFSPINOR" statement
 *
 * - the numbers represent the certain implementation of Hopping_Matrix
 *
 * #if defined _USE_HALFSPINOR
 *   #if ((defined SSE2)||(defined SSE3))
 *     1. 
 *    #elif (defined BGL && defined XLC)
 *     2.
 *   #else
 *     3.
 *   #endif
 * #else * thats _USE_HALFSPINOR *
 *   #if ((defined SSE2)||(defined SSE3))
 *     #if (defined _USE_TSPLITPAR)
 *     4.
 *     #else
 *     5.
 *     #endif
 *   #elif (defined BGL && defined XLC)
 *     6.
 *   #elif defined XLC
 *     7.
 *   * else of If defined SSE2  and if defined XLC *
 *   #else
 *     8.
 *   #endif
 * #endif * thats _USE_HALFSPINOR *
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
#include "global.h"
#include "su3.h"
#ifdef MPI
#  include "xchange_field.h"
#  ifdef _USE_TSPLITPAR
#    include "xchange_field_tslice.h"
#  endif
#  if defined _USE_HALFSPINOR
#    include "xchange_halffield.h"
#  endif
#endif
#include "boundary.h"
#include "init_dirac_halfspinor.h"
#include "update_backward_gauge.h"
#include "Hopping_Matrix.h"

#if defined _USE_HALFSPINOR
#  if ((defined SSE2)||(defined SSE3))
#    include "sse.h"
#    include "operator/halfspinor_sse_dbl.c"

#  elif (defined BGL && defined XLC)
#    include "bgl.h"
#    include "operator/halfspinor_bg_dbl.c"

#  elif (defined BGQ && defined XLC)
#    include "bgq.h"
#    include "operator/halfspinor_bgq_dbl.c"

#  else
#    include "operator/halfspinor_dbl.c"

#  endif

#else /* thats _USE_HALFSPINOR */

#  if (((defined SSE2)||(defined SSE3)) && defined _USE_TSPLITPAR)
#    include "sse.h"
#    include "operator/hopping_sse_dbl.c"

#  elif (defined BGL && defined XLC)
#    include "bgl.h"
#    include "operator/hopping_bg_dbl.c"

#  else
#    if ((defined SSE2)||(defined SSE3))
#      include "sse.h"

#    elif (defined BGQ && defined XLC)
#      include "bgq.h"
#      include "xlc_prefetch.h"

#    elif defined XLC
#      include"xlc_prefetch.h"

#    endif
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {

#    include "operator/hopping_body_dbl.c"

  return;
}
#  endif

#endif /* thats _USE_HALFSPINOR */

