/***********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "monomial/monomial.h"
#include "operator/Hopping_Matrix.h"
#include "gettime.h"
#include "clovernd_trlog_monomial.h"

void clovernd_trlog_derivative(const int id, hamiltonian_field_t * const hf) {
  //monomial * mnl = &monomial_list[id];
  /* this term has no derivative */
  /* so a dummy function         */
  if(g_proc_id == 0 && g_debug_level > 4) {
    printf("called clovernd_trlog_derivative for id %d, which is a dummy function\n", id);
  }
  return;
}


void clovernd_trlog_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  mnl->energy0 = 0.;

  init_sw_fields();
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  /*compute the contribution from the clover trlog term */
  mnl->energy0 = -sw_trace_nd(EE, mnl->mubar, mnl->epsbar);
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called clovernd_trlog_heatbath for id %d E = %e\n", id, mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}

double clovernd_trlog_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  mnl->energy1 = 0.;
  sw_term( (const su3**) hf->gaugefield, mnl->kappa, mnl->c_sw); 
  /*compute the contribution from the clover trlog term */
  mnl->energy1 = -sw_trace_nd(EE, mnl->mubar, mnl->epsbar);
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called clovernd_trlog_acc for id %d dH = %1.10e\n", 
	   id, mnl->energy1 - mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return(mnl->energy1 - mnl->energy0);
}
