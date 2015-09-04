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
#include "gettime.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "measure_gauge_action.h"
#include "measure_rectangles.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "gauge_monomial.h"

/* this function calculates the derivative of the momenta: equation 13 of Gottlieb */
void gauge_derivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double factor = -1. * g_beta/3.0;
  if(mnl->use_rectangles) {
    mnl->forcefactor = 1.;
    factor = -mnl->c0 * g_beta/3.0;
  }
  
  double atime, etime;
  atime = gettime();
#ifdef OMP
#pragma omp parallel
  {
#endif

  su3 ALIGN v, w;
  int i, mu;
  su3 *z;
  su3adj *xm;

#ifdef OMP
#pragma omp for
#endif
  for(i = 0; i < VOLUME; i++) { 
    for(mu=0;mu<4;mu++) {
      z=&hf->gaugefield[i][mu];
      xm=&hf->derivative[i][mu];
      get_staples(&v,i,mu, (const su3**) hf->gaugefield); 
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
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}

/* this function calculates the derivative of the momenta: equation 13 of Gottlieb */
void gauge_EMderivative(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double factor = -1. * g_beta/3.0;
  if(mnl->use_rectangles) {
    mnl->forcefactor = 1.;
    factor = -mnl->c0 * g_beta/3.0;
  }
  
  double atime, etime;
  atime = gettime();
#ifdef OMP
#pragma omp parallel
  {
#endif

  su3 ALIGN v, w;
  int i, mu;
  su3 *z;
  su3adj *xm;

#ifdef OMP
#pragma omp for
#endif
  for(i = 0; i < VOLUME; i++) { 
    // electric part
    z=&hf->gaugefield[i][0];
    xm=&hf->derivative[i][0];
    get_staples(&v, i, 0, (const su3**) hf->gaugefield); 
    _su3_times_su3d(w,*z,v);
    _trace_lambda_mul_add_assign((*xm), (1.+mnl->glambda)*factor, w);
    // lambda only acts on the plaquette, effectively changing c0 in the spatial and temporal parts, c1 remains untouched
    if(mnl->use_rectangles) {
	    get_rectangle_staples(&v, i, 0);
	    _su3_times_su3d(w, *z, v);
	    _trace_lambda_mul_add_assign((*xm), factor*mnl->c1/mnl->c0, w);
    }
    // magnetic part
    for(mu=1;mu<4;mu++) {
      z=&hf->gaugefield[i][mu];
      xm=&hf->derivative[i][mu];

      get_spacelike_staples(&v, i, mu, (const su3**) hf->gaugefield); 
      _su3_times_su3d(w, *z, v);
      _trace_lambda_mul_add_assign((*xm), (1.-mnl->glambda)*factor, w);

      get_timelike_staples(&v, i, mu, (const su3**) hf->gaugefield); 
      _su3_times_su3d(w, *z, v);
      _trace_lambda_mul_add_assign((*xm), (1.+mnl->glambda)*factor, w);
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
  etime = gettime();
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# Time for %s monomial derivative: %e s\n", mnl->name, etime-atime);
  }
  return;
}

void gauge_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();
  if(mnl->use_rectangles) mnl->c0 = 1. - 8.*mnl->c1;
  
  mnl->energy0 = g_beta*(mnl->c0 * measure_gauge_action( (const su3**) hf->gaugefield, mnl->glambda));
  if(mnl->use_rectangles) {
    mnl->energy0 += g_beta*(mnl->c1 * measure_rectangles( (const su3**) hf->gaugefield));
  }
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial heatbath: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called gauge_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  return;
}

double gauge_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  double atime, etime;
  atime = gettime();
  mnl->energy1 = g_beta*(mnl->c0 * measure_gauge_action( (const su3**) hf->gaugefield, mnl->glambda));
  if(mnl->use_rectangles) {
    mnl->energy1 += g_beta*(mnl->c1 * measure_rectangles( (const su3**) hf->gaugefield));
  }
  etime = gettime();
  if(g_proc_id == 0) {
    if(g_debug_level > 1) {
      printf("# Time for %s monomial acc step: %e s\n", mnl->name, etime-atime);
    }
    if(g_debug_level > 3) {
      printf("called gauge_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy0 - mnl->energy1);
    }
  }
  return(mnl->energy0 - mnl->energy1);
}
