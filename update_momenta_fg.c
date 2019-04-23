/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *
 * some changes by C. Urbach 2002-2008,2012
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <lime.h>
#include "sighandler.h"
#include "read_input.h"
#include "monomial/monomial.h"
#include "init/init_gauge_fg.h"
#include "operator/clover_leaf.h"

#include "global.h"
#include "gettime.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "expo.h"
#include "sse.h"
#include "xchange/xchange.h"
#include "hamiltonian_field.h"
#include "init/init_gauge_field.h"
#ifdef DDalphaAMG
#include "DDalphaAMG_interface.h"
#endif
/*******************************************************
 *
 * Temporarily updates the gauge field corresponding to 
 * the approximated force gradient term to finally update
 * the momenta
 *
 *******************************************************/

void update_momenta_fg(int * mnllist, double step, const int no,
		       hamiltonian_field_t * const hf, double step0) {
  double atime, etime;
  atime = gettime();
#ifdef DDalphaAMG
  MG_update_gauge(0.0);
#endif

  /* #ifdef TM_USE_OMP
     #define static
     #pragma omp parallel
     {
     #endif
  */

  int i,mu;
  double step_fg;
  static su3 v,w;
  su3 *z;
  su3 *ztmp;
  static su3adj deriv;
  su3adj *Fm;

  step_fg=-step0*step0/24;
  /*
     #ifdef _KOJAK_INST
     #pragma pomp inst begin(updategauge)
     #endif

     #ifdef TM_USE_OMP
     #pragma omp parallel for
     #endif
  */

  for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand);i++) {
    for(int mu=0;mu<4;mu++) {
      _zero_su3adj(hf->derivative[i][mu]);
    }
  }

  for(int k = 0; k < no; k++) {
    if(monomial_list[ mnllist[k] ].derivativefunction != NULL) {
      monomial_list[ mnllist[k] ].derivativefunction(mnllist[k], hf);
    }
  }

#ifdef TM_USE_MPI
  xchange_deri(hf->derivative);
#endif


  /* #ifdef TM_USE_OMP
     #pragma omp parallel for
     #endif
  */

  for(i = 0; i < VOLUME; i++) { 
    for(mu = 0; mu < 4; mu++){
      /* Cope gauge field to be temporarily updated */
      z = &hf->gaugefield[i][mu];
      ztmp = &gauge_fg[i][mu];
      _su3_assign(*ztmp,*z);  
 
      /* Calculate approximated force gradient term and update temporary gauge field */
      Fm = &hf->derivative[i][mu];
      _zero_su3adj(deriv);
      _su3adj_assign_const_times_su3adj(deriv, step_fg, *Fm);
      /*_su3adj_assign_const_times_su3adj(deriv, 0.0, *Fm);*/
      exposu3(&w,&deriv);
      restoresu3(&v,&w);
      _su3_times_su3(w, v, *z);
      _su3_assign(*z, w);
    }
  }

#ifdef TM_USE_MPI
     /* for parallelization */
     xchange_gauge(hf->gaugefield);
#endif
#ifdef DDalphaAMG
     MG_update_gauge(0.0);
#endif

   /*Convert to a 32 bit gauge field, after xchange*/
   convert_32_gauge_field(g_gauge_field_32, hf->gaugefield, VOLUMEPLUSRAND + g_dbw2rand);

   /* The backward copy of gaugefield is not updated here! */
   hf->update_gauge_copy = 1;
   g_update_gauge_copy = 1;
   g_update_gauge_copy_32 = 1;


   /* #ifdef TM_USE_OMP
      #pragma omp parallel for
      #endif
   */
   /* Calculate derivate based on the temporary updated
      gauge field U'=ztmp:
      1) Set derivative to zero
      2) Recalcuate derivate
   */
    
  for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand);i++) {
    for(int mu=0;mu<4;mu++) {
      _zero_su3adj(hf->derivative[i][mu]);
    }
  }

  for(int k = 0; k < no; k++) {
    if(monomial_list[ mnllist[k] ].derivativefunction != NULL) {
      monomial_list[ mnllist[k] ].derivativefunction(mnllist[k], hf);
    }
  }

#ifdef TM_USE_MPI
  xchange_deri(hf->derivative);
#endif

  for(i = 0; i < VOLUME; i++) { 
    for(mu = 0; mu < 4; mu++){
      /* Update momenta (the minus comes from an extra minus in trace_lambda)
	 and restore initial gauge field */
      _su3adj_minus_const_times_su3adj(hf->momenta[i][mu], step, hf->derivative[i][mu]);

      z = &hf->gaugefield[i][mu];
      ztmp = &gauge_fg[i][mu];
      _su3_assign(*z,*ztmp);

    }
  }

  /* #ifdef TM_USE_OMP
     } /* OpenMP parallel closing brace /
     #endif
  */
  
#ifdef TM_USE_MPI
  /* for parallelization */
  xchange_gauge(hf->gaugefield);
#endif
#ifdef DDalphaAMG
  MG_update_gauge(0.0);
#endif

  /*Convert to a 32 bit gauge field, after xchange*/
  convert_32_gauge_field(g_gauge_field_32, hf->gaugefield, VOLUMEPLUSRAND + g_dbw2rand);
  
  /*
   * The backward copy of the gauge field
   * is not updated here!
   */
  hf->update_gauge_copy = 1;
  g_update_gauge_copy = 1;
  g_update_gauge_copy_32 = 1;


  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time gauge update: %e s\n", etime-atime); 
  } 
  return;

  /* #ifdef _KOJAK_INST
     #pragma pomp inst end(updategauge)
     #endif
  */
}
