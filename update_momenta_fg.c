/***********************************************************************
 *
 * Copyright (C) 2017 Jacob Finkenrath
 *               2018 Bartosz Kostrzewa
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
# include<config.h>
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

inline void calculate_fg(const double step_fg,
                         hamiltonian_field_t * const hf){
#ifdef TM_USE_OMP
#define static
#pragma omp parallel
  {
#endif

  static su3 v,w;
  su3 *z;
  su3 *ztmp;
  static su3adj deriv;
  su3adj *Fm;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(int i = 0; i < VOLUME; i++) { 
    for(int mu = 0; mu < 4; mu++){
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
      restoresu3(&v,&w);
      _su3_assign(*z, v);
    }
  }
#ifdef TM_USE_OMP
  } // OpenMP parallel section closing brace
#undef static
#endif
}

inline void fg_update_momenta_reset_gaugefield(const double step,
                                               hamiltonian_field_t * const hf){
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
  su3 *z;
  su3 *ztmp;
#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(int i = 0; i < VOLUME; i++) { 
    for(int mu = 0; mu < 4; mu++){
      /* Update momenta (the minus comes from an extra minus in trace_lambda)
       and restore initial gauge field */
      _su3adj_minus_const_times_su3adj(hf->momenta[i][mu], step, hf->derivative[i][mu]);
  
      z = &hf->gaugefield[i][mu];
      ztmp = &gauge_fg[i][mu];
      _su3_assign(*z,*ztmp);
  
    }
  }
#ifdef TM_USE_OMP
  } // OpenMP parallel section closing brace
#endif
}

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
  if (g_exposu3_no_c == 0) init_exposu3();

  double step_fg=-step0*step0/24;

#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand);i++) {
    for(int mu=0;mu<4;mu++) {
      _zero_su3adj(hf->derivative[i][mu]);
    }
  }

  // calculate derivatives to estimate force gradient
  for(int k = 0; k < no; k++) {
    if(monomial_list[ mnllist[k] ].derivativefunction != NULL) {
      monomial_list[ mnllist[k] ].derivativefunction(mnllist[k], hf);
    }
  }

#ifdef TM_USE_MPI
  xchange_deri(hf->derivative);
#endif
  // estimate force gradient and propagate to gauge field
  calculate_fg(step_fg, hf);

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

  // calculate forces with force-gradient updated gauge field
#ifdef TM_USE_OMP
#pragma omp parallel for
#endif
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
  
  // and finally update the momenta and reset the gauge field 
  fg_update_momenta_reset_gaugefield(step, hf);

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
}
