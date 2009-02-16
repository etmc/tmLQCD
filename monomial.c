/***********************************************************************
 * $Id$ 
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
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "ranlxd.h"
#include "sse.h"
#include "linalg_eo.h"
#include "default_input_values.h"
#include "read_input.h"
#include "monomial.h"



monomial monomial_list[max_no_monomials];
int no_monomials = 0;
int no_gauge_monomials = 0;
int no_ndpoly_monomials = 0;
static spinor * _pf;

int add_monomial(const int type) {
  
  if(no_monomials == max_no_monomials) {
    fprintf(stderr, "maximal number of monomials %d exceeded!\n", max_no_monomials);
    exit(-1);
  }
  if(type == DET) {
    monomial_list[no_monomials].hbfunction = &det_heatbath;
    monomial_list[no_monomials].accfunction = &det_acc;
    monomial_list[no_monomials].derivativefunction = &det_derivative;
  }
  else if(type == DETRATIO) {
    monomial_list[no_monomials].hbfunction = &detratio_heatbath;
    monomial_list[no_monomials].accfunction = &detratio_acc;
    monomial_list[no_monomials].derivativefunction = &detratio_derivative;
  }
  else if(type == GAUGE) {
    if(no_gauge_monomials > 0) {
      fprintf(stderr, "maximal number of gauge monomials exceeded!\n");
      exit(-1);
    }
    monomial_list[no_monomials].hbfunction = &gauge_heatbath;
    monomial_list[no_monomials].accfunction = &gauge_acc;
    monomial_list[no_monomials].derivativefunction = &gauge_derivative;
    no_gauge_monomials++;
    if( monomial_list[no_monomials].gtype == 0) {
      monomial_list[no_monomials].use_rectangles = 0;
      monomial_list[no_monomials].c0 = 1.;
      monomial_list[no_monomials].c1 = 0.;
    }
    if(monomial_list[no_monomials].gtype == 1) {
      monomial_list[no_monomials].use_rectangles = 1;
      monomial_list[no_monomials].c1 = -0.083333333;
      monomial_list[no_monomials].c0 = 1.-8*monomial_list[no_monomials].c1;
    }
    if(monomial_list[no_monomials].gtype == 2) {
      monomial_list[no_monomials].use_rectangles = 1;
      monomial_list[no_monomials].c1 = -0.331;
      monomial_list[no_monomials].c0 = 1.-8*monomial_list[no_monomials].c1;
    }
    if(monomial_list[no_monomials].gtype == 4) {
      monomial_list[no_monomials].use_rectangles = 1;
      monomial_list[no_monomials].c1 = -1.4088;
      monomial_list[no_monomials].c0 = 1.-8*monomial_list[no_monomials].c1;
    }
    if(!monomial_list[no_monomials].use_rectangles) {
      g_rgi_C1 = 0.;
    }
    else g_rgi_C1 = monomial_list[no_monomials].c1;
  }
  else if(type == NDPOLY) {
    if(no_ndpoly_monomials > 0) {
      fprintf(stderr, "maximal number of ndpoly monomials (1) exceeded!\n");
      exit(-1);
    }
    monomial_list[no_monomials].hbfunction = &ndpoly_heatbath;
    monomial_list[no_monomials].accfunction = &ndpoly_acc;
    monomial_list[no_monomials].derivativefunction = &ndpoly_derivative;
    no_ndpoly_monomials++;
  }
  else {
    fprintf(stderr, "Unknown monomial type!\n");
    return(-1);
  }
  monomial_list[no_monomials].pf = NULL;
  monomial_list[no_monomials].csg_field = NULL;
  monomial_list[no_monomials].csg_field2 = NULL;
  monomial_list[no_monomials].csg_index_array = NULL;
  monomial_list[no_monomials].csg_index_array2 = NULL;
  monomial_list[no_monomials].csg_N = 0;
  monomial_list[no_monomials].csg_N2 = 0;
  monomial_list[no_monomials].csg_n = 1;
  monomial_list[no_monomials].csg_n2 = 1;
  monomial_list[no_monomials].kappa = _default_g_kappa;
  monomial_list[no_monomials].kappa2 = _default_g_kappa;
  monomial_list[no_monomials].mu = _default_g_mu;
  monomial_list[no_monomials].mu2 = _default_g_mu;
  monomial_list[no_monomials].epsilon = _default_g_epsbar;
  monomial_list[no_monomials].timescale = _default_timescale;
  monomial_list[no_monomials].accprec = _default_g_eps_sq_acc;
  monomial_list[no_monomials].forceprec = _default_g_eps_sq_force;
  monomial_list[no_monomials].maxiter = _default_max_solver_iterations;
  monomial_list[no_monomials].solver = _default_solver_flag;
  monomial_list[no_monomials].even_odd_flag = _default_even_odd_flag;
  monomial_list[no_monomials].forcefactor = 1.;
  monomial_list[no_monomials].use_rectangles = 0;
  monomial_list[no_monomials].c1 = _default_g_rgi_C1;
  monomial_list[no_monomials].c0 = 1.;
  monomial_list[no_monomials].beta = _default_g_beta;  
  monomial_list[no_monomials].rngrepro = _default_reproduce_randomnumber_flag;
  monomial_list[no_monomials].initialised = 1;

  no_monomials++;
  return(no_monomials);
}


int init_monomials(const int V, const int even_odd_flag) {
  int i, no=0;
  spinor * __pf = NULL;
  for(i = 0; i < no_monomials; i++) {
    if(monomial_list[i].type != GAUGE) no++;
  }
  if(no_monomials > 0) {
    if((void*)(_pf = (spinor*)calloc(no*V+1, sizeof(spinor))) == NULL) {
      printf ("malloc errno in monomial pf fields: %d\n",errno); 
      errno = 0;
      return(1);
    }
    else {
#if ( defined SSE || defined SSE2 || defined SSE3)
      __pf = (spinor*)(((unsigned long int)(_pf)+ALIGN_BASE)&~ALIGN_BASE);
#else
      __pf = _pf;
#endif
    }
  }

  no = 0;
  for(i = 0; i < no_monomials; i++) {
    if(monomial_list[i].type != GAUGE) {
      monomial_list[i].pf = __pf+no*V;
      no++;
      monomial_list[i].rngrepro = reproduce_randomnumber_flag;
    }
    else {
      monomial_list[i].pf = NULL;
      if(!monomial_list[i].use_rectangles) {
	monomial_list[i].c1 = 0.;
	g_rgi_C1 = 0.;
	monomial_list[i].c0 = 1.;
	g_rgi_C0 = 1.;
      }
      else {
	monomial_list[i].c0 = 1. - 8.*monomial_list[i].c1;
      }
    }
    monomial_list[i].id = i;
    monomial_list[i].even_odd_flag = even_odd_flag;
  }
  return(0);
}

void free_monomials() {
  
  free(_pf);
  return;
}

void dummy_heatbath(const int id) {
  return;
}

double dummy_acc(const int id) {
  return(0.);
}
