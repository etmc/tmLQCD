/***********************************************************************
 *
 * Copyright (C) 2007,2008 Jan Volkholz, Carsten Urbach
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
#include <tmlqcd_config.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "boundary.h"
#include "deriv_Sb_D_psi.h"
#include "gettime.h"
#include "global.h"
#include "hamiltonian_field.h"
#include "su3.h"
#include "xchange/xchange.h"

/*----------------------------------------------------------------------------*/

void deriv_Sb_D_psi(spinor* const l, spinor* const k, hamiltonian_field_t* const hf,
                    const double factor) {
  tm_stopwatch_push(&g_timers, __func__, "");
  /* for parallelization */
#ifdef TM_USE_MPI
  xchange_lexicfield(k);
  xchange_lexicfield(l);
#endif

#ifdef TM_USE_OMP
#define static
#pragma omp parallel
  {
#endif

    int ix, iy;
    su3* restrict up ALIGN;
    su3* restrict um ALIGN;
    static su3 v1, v2;
    static su3_vector psia, psib, phia, phib;
    static spinor rr;
    /*   spinor * restrict r ALIGN; */
    spinor* restrict sp ALIGN;
    spinor* restrict sm ALIGN;

#ifdef TM_USE_OMP
#undef static
#endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(derivSb)
#endif

    /************** loop over all lattice sites ****************/
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for (ix = 0; ix < (VOLUME); ix++) {
      rr = (*(l + ix));
      /*     rr=g_spinor_field[l][icx-ioff]; */

      /*multiply the left vector with gamma5*/
      _vector_minus_assign(rr.s2, rr.s2);
      _vector_minus_assign(rr.s3, rr.s3);

      /*********************** direction +0 ********************/

      iy = g_iup[ix][0];

      sp = k + iy;
      up = &hf->gaugefield[ix][0];

      _vector_add(psia, (*sp).s0, (*sp).s2);
      _vector_add(psib, (*sp).s1, (*sp).s3);

      _vector_add(phia, rr.s0, rr.s2);
      _vector_add(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka0, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][0], 2. * factor, v1);

      /************** direction -0 ****************************/

      iy = g_idn[ix][0];

      sm = k + iy;
      um = &hf->gaugefield[iy][0];

      _vector_sub(psia, (*sm).s0, (*sm).s2);
      _vector_sub(psib, (*sm).s1, (*sm).s3);

      _vector_sub(phia, rr.s0, rr.s2);
      _vector_sub(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka0, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][0], 2. * factor, v1);

      /*************** direction +1 **************************/

      iy = g_iup[ix][1];

      sp = k + iy;
      up = &hf->gaugefield[ix][1];

      _vector_i_add(psia, (*sp).s0, (*sp).s3);
      _vector_i_add(psib, (*sp).s1, (*sp).s2);

      _vector_i_add(phia, rr.s0, rr.s3);
      _vector_i_add(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka1, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][1], 2. * factor, v1);

      /**************** direction -1 *************************/

      iy = g_idn[ix][1];

      sm = k + iy;
      um = &hf->gaugefield[iy][1];

      _vector_i_sub(psia, (*sm).s0, (*sm).s3);
      _vector_i_sub(psib, (*sm).s1, (*sm).s2);

      _vector_i_sub(phia, rr.s0, rr.s3);
      _vector_i_sub(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka1, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][1], 2. * factor, v1);

      /*************** direction +2 **************************/

      iy = g_iup[ix][2];

      sp = k + iy;
      up = &hf->gaugefield[ix][2];

      _vector_add(psia, (*sp).s0, (*sp).s3);
      _vector_sub(psib, (*sp).s1, (*sp).s2);

      _vector_add(phia, rr.s0, rr.s3);
      _vector_sub(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka2, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][2], 2. * factor, v1);

      /***************** direction -2 ************************/

      iy = g_idn[ix][2];

      sm = k + iy;
      um = &hf->gaugefield[iy][2];

      _vector_sub(psia, (*sm).s0, (*sm).s3);
      _vector_add(psib, (*sm).s1, (*sm).s2);

      _vector_sub(phia, rr.s0, rr.s3);
      _vector_add(phib, rr.s1, rr.s2);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka2, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][2], 2. * factor, v1);

      /****************** direction +3 ***********************/

      iy = g_iup[ix][3];

      sp = k + iy;
      up = &hf->gaugefield[ix][3];

      _vector_i_add(psia, (*sp).s0, (*sp).s2);
      _vector_i_sub(psib, (*sp).s1, (*sp).s3);

      _vector_i_add(phia, rr.s0, rr.s2);
      _vector_i_sub(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, phia, psia, phib, psib);
      _su3_times_su3d(v2, *up, v1);
      _complex_times_su3(v1, ka3, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[ix][3], 2. * factor, v1);

      /***************** direction -3 ************************/

      iy = g_idn[ix][3];

      sm = k + iy;
      um = &hf->gaugefield[iy][3];

      _vector_i_sub(psia, (*sm).s0, (*sm).s2);
      _vector_i_add(psib, (*sm).s1, (*sm).s3);

      _vector_i_sub(phia, rr.s0, rr.s2);
      _vector_i_add(phib, rr.s1, rr.s3);

      _vector_tensor_vector_add(v1, psia, phia, psib, phib);
      _su3_times_su3d(v2, *um, v1);
      _complex_times_su3(v1, ka3, v2);
      _trace_lambda_mul_add_assign_nonlocal(hf->derivative[iy][3], 2. * factor, v1);

      /****************** end of loop ************************/
    }
#ifdef _KOJAK_INST
#pragma pomp inst end(derivSb)
#endif

#ifdef TM_USE_OMP
  } /*OpenMP closing brace */
#endif
  tm_stopwatch_pop(&g_timers, 0, 1, "");
}
