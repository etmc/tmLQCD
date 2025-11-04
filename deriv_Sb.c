/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
 *
 * some changes to initial version by Carsten Urbach
 *
 * BG version Copyright (C) 2006, 2007 Carsten Urbach
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
 * deriv_Sb: function to compute the derivative
 * of the phi^{\dag} Q psi with respect
 * to the generators of the gauge group.
 * without the clover part.
 *
 * Author: Martin Hasenbusch <Martin.Hasenbusch@desy.de>
 * Date: Fri Oct 26 15:06:27 MEST 2001
 *
 *  both l and k are input
 *  for ieo = 0
 *  l resides on even lattice points and k on odd lattice points
 *  for ieo = 1
 *  l resides on odd lattice points and k on even lattice points
 *  the output is a su3adj field that is written to df0[][]
 *
 ************************************************************************/

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "boundary.h"
#include "deriv_Sb.h"
#include "gettime.h"
#include "global.h"
#include "hamiltonian_field.h"
#include "su3.h"
#include "update_backward_gauge.h"
#include "xchange/xchange.h"

void deriv_Sb(const int ieo, spinor* const l, spinor* const k, hamiltonian_field_t* const hf,
              const double factor) {
  tm_stopwatch_push(&g_timers, __func__, "");
#ifdef _GAUGE_COPY
  if (g_update_gauge_copy) {
    update_backward_gauge(hf->gaugefield);
  }
#endif
  /* for parallelization */
#ifdef TM_USE_MPI
  xchange_2fields(k, l, ieo);
#endif

#ifdef TM_USE_OMP
#define static
#pragma omp parallel
  {
#endif
    int ix, iy;
    int ioff, icx, icy;
    su3* restrict up ALIGN;
    su3* restrict um ALIGN;
    static su3 v1, v2;
    static su3_vector psia, psib, phia, phib;
    static spinor rr;
    spinor* restrict sp ALIGN;
    spinor* restrict sm ALIGN;

#ifdef TM_USE_OMP
#undef static
#endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(derivSb)
#endif

#ifdef BGL
    __alignx(16, l);
    __alignx(16, k);
#endif

    if (ieo == 0) {
      ioff = 0;
    } else {
      ioff = (VOLUME + RAND) / 2;
    }

    /************** loop over all lattice sites ****************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (icx = ioff; icx < (VOLUME / 2 + ioff); icx++) {
      ix = g_eo2lexic[icx];
      rr = (*(l + (icx - ioff)));
      /*     rr=g_spinor_field[l][icx-ioff]; */

      /*multiply the left vector with gamma5*/
      _vector_minus_assign(rr.s2, rr.s2);
      _vector_minus_assign(rr.s3, rr.s3);

      /*********************** direction +0 ********************/

      iy = g_iup[ix][0];
      icy = g_lexic2eosub[iy];

      sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      up = &g_gauge_field_copy[icx][0];
#else
    up = &hf->gaugefield[ix][0];
#endif
      _vector_add(psia, sp->s0, sp->s2);
      _vector_add(psib, sp->s1, sp->s3);

      _vector_add(phia, rr.s0, rr.s2);
      _vector_add(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka0, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][0], 2. * factor, v1);

      /************** direction -0 ****************************/

      iy = g_idn[ix][0];
      icy = g_lexic2eosub[iy];

      sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      um = up + 1;
#else
    um = &hf->gaugefield[iy][0];
#endif

      _vector_sub(psia, sm->s0, sm->s2);
      _vector_sub(psib, sm->s1, sm->s3);

      _vector_sub(phia, rr.s0, rr.s2);
      _vector_sub(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka0, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][0], 2. * factor, v1);

      /*************** direction +1 **************************/

      iy = g_iup[ix][1];
      icy = g_lexic2eosub[iy];

      sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      up = um + 1;
#else
    up = &hf->gaugefield[ix][1];
#endif
      _vector_i_add(psia, sp->s0, sp->s3);
      _vector_i_add(psib, sp->s1, sp->s2);

      _vector_i_add(phia, rr.s0, rr.s3);
      _vector_i_add(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka1, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][1], 2. * factor, v1);

      /**************** direction -1 *************************/

      iy = g_idn[ix][1];
      icy = g_lexic2eosub[iy];

      sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      um = up + 1;
#else
    um = &hf->gaugefield[iy][1];
#endif
      _vector_i_sub(psia, sm->s0, sm->s3);
      _vector_i_sub(psib, sm->s1, sm->s2);

      _vector_i_sub(phia, rr.s0, rr.s3);
      _vector_i_sub(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka1, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][1], 2. * factor, v1);

      /*************** direction +2 **************************/

      iy = g_iup[ix][2];
      icy = g_lexic2eosub[iy];

      sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      up = um + 1;
#else
    up = &hf->gaugefield[ix][2];
#endif
      _vector_add(psia, sp->s0, sp->s3);
      _vector_sub(psib, sp->s1, sp->s2);

      _vector_add(phia, rr.s0, rr.s3);
      _vector_sub(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka2, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][2], 2. * factor, v1);

      /***************** direction -2 ************************/

      iy = g_idn[ix][2];
      icy = g_lexic2eosub[iy];

      sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      um = up + 1;
#else
    um = &hf->gaugefield[iy][2];
#endif
      _vector_sub(psia, sm->s0, sm->s3);
      _vector_add(psib, sm->s1, sm->s2);

      _vector_sub(phia, rr.s0, rr.s3);
      _vector_add(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka2, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][2], 2. * factor, v1);

      /****************** direction +3 ***********************/

      iy = g_iup[ix][3];
      icy = g_lexic2eosub[iy];

      sp = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      up = um + 1;
#else
    up = &hf->gaugefield[ix][3];
#endif
      _vector_i_add(psia, sp->s0, sp->s2);
      _vector_i_sub(psib, sp->s1, sp->s3);

      _vector_i_add(phia, rr.s0, rr.s2);
      _vector_i_sub(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka3, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][3], 2. * factor, v1);

      /***************** direction -3 ************************/

      iy = g_idn[ix][3];
      icy = g_lexic2eosub[iy];

      sm = k + icy;
#if (defined _GAUGE_COPY && !defined _USE_HALFSPINOR && !defined _USE_TSPLITPAR)
      um = up + 1;
#else
    um = &hf->gaugefield[iy][3];
#endif
      _vector_i_sub(psia, sm->s0, sm->s2);
      _vector_i_add(psib, sm->s1, sm->s3);

      _vector_i_sub(phia, rr.s0, rr.s2);
      _vector_i_add(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka3, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][3], 2. * factor, v1);

      /****************** end of loop ************************/
    }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch_pop(&g_timers, 0, 1, "");
#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif
}
