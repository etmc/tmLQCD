/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
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
#include <stdio.h>
#include <math.h>
#ifdef OMP 
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "sse.h"
#include "start.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "measure_gauge_action.h"
#include "measure_rectangles.h"
#include "monomial.h"
#include "hamiltonian_field.h"
#include "gauge_monomial.h"

/* this function calculates the derivative of the momenta: equation 13 of Gottlieb */
void gauge_derivative(const int id, hamiltonian_field_t * const hf) {
#ifdef OMP
#define static
#pragma omp parallel
  {
#endif

  static su3 v, w;
  int i, mu;
  su3 *z;
  su3adj *xm;
  monomial * mnl = &monomial_list[id];
  double factor = -1. * g_beta/3.0;

#ifdef OMP
#undef static
#endif

  if(mnl->use_rectangles) {
    mnl->forcefactor = 1.;
    factor = -mnl->c0 * g_beta/3.0;
  }
  
#ifdef OMP
#pragma omp for
#endif
  for(i = 0; i < VOLUME; i++) { 
    for(mu=0;mu<4;mu++) {
      z=&hf->gaugefield[i][mu];
      xm=&hf->derivative[i][mu];
      v=get_staples(i,mu, hf->gaugefield); 
      _su3_times_su3d(w,*z,v);
      _trace_lambda_mul_add_assign((*xm), factor, w);
      
      if(mnl->use_rectangles) {
	get_rectangle_staples(&v, i, mu);
	_su3_times_su3d(w, *z, v);
	_trace_lambda_mul_add_assign((*xm), factor*mnl->c1/mnl->c0, w);
      }
    }
  }

#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}

void gauge_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  
  if(mnl->use_rectangles) mnl->c0 = 1. - 8.*mnl->c1;
  
  mnl->energy0 = g_beta*(mnl->c0 * measure_gauge_action(hf->gaugefield));
  if(mnl->use_rectangles) {
    mnl->energy0 += g_beta*(mnl->c1 * measure_rectangles(hf->gaugefield));
  }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called gauge_heatbath for id %d %d\n", id, mnl->even_odd_flag);
  }
}

double gauge_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  
  mnl->energy1 = g_beta*(mnl->c0 * measure_gauge_action(hf->gaugefield));
  if(mnl->use_rectangles) {
    mnl->energy1 += g_beta*(mnl->c1 * measure_rectangles(hf->gaugefield));
    }
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called gauge_acc for id %d %d dH = %1.4e\n", 
	   id, mnl->even_odd_flag, mnl->energy0 - mnl->energy1);
  }
  return(mnl->energy0 - mnl->energy1);
}
