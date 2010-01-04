/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Carsten Urbach
 *
 * Modified by Jenifer Gonzalez Lopez 2009/03/31
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
  monomial_list[no_monomials].hbfunction = &dummy_heatbath;
  monomial_list[no_monomials].accfunction = &dummy_acc;
  monomial_list[no_monomials].derivativefunction = &dummy_derivative;

  monomial_list[no_monomials].pf = NULL;
  monomial_list[no_monomials].pf2 = NULL;
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
  monomial_list[no_monomials].mubar = _default_g_mubar;
  monomial_list[no_monomials].mubar2 = _default_g_mubar;
  monomial_list[no_monomials].epsbar = _default_g_epsbar;
  monomial_list[no_monomials].epsbar2 = _default_g_epsbar;
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
  monomial_list[no_monomials].eta = _default_g_eta;
  monomial_list[no_monomials].ct = _default_g_Ct; 
  monomial_list[no_monomials].cs = _default_g_Cs;
  monomial_list[no_monomials].c1ss = _default_g_C1ss; 
  monomial_list[no_monomials].c1tss = _default_g_C1tss; 
  monomial_list[no_monomials].c1tts = _default_g_C1tts; 
  monomial_list[no_monomials].rngrepro = _default_reproduce_randomnumber_flag;
  monomial_list[no_monomials].initialised = 1;
  if(monomial_list[no_monomials].type == NDDETRATIO) {
    monomial_list[no_monomials].timescale = -5;
  }

  no_monomials++;
  return(no_monomials);
}


int init_monomials(const int V, const int even_odd_flag) {
  int i, no=0;
  spinor * __pf = NULL;
  for(i = 0; i < no_monomials; i++) {
    if((monomial_list[i].type != GAUGE) && (monomial_list[i].type != SFGAUGE)) no++;
    /* non-degenerate monomials need two pseudo fermion fields */
    if((monomial_list[i].type == NDPOLY) || (monomial_list[i].type == NDDETRATIO)) no++;
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
    
    if((monomial_list[i].type != GAUGE) && (monomial_list[i].type != SFGAUGE)) {
          
      monomial_list[i].pf = __pf+no*V;
      no++;
      monomial_list[i].rngrepro = reproduce_randomnumber_flag;

      if(monomial_list[i].type == DET) {
	monomial_list[i].hbfunction = &det_heatbath;
	monomial_list[i].accfunction = &det_acc;
	monomial_list[i].derivativefunction = &det_derivative;
      }
      else if(monomial_list[i].type == DETRATIO) {
	monomial_list[i].hbfunction = &detratio_heatbath;
	monomial_list[i].accfunction = &detratio_acc;
	monomial_list[i].derivativefunction = &detratio_derivative;
      }
      else if(monomial_list[i].type == NDPOLY) {
	if(no_ndpoly_monomials > 0) {
	  fprintf(stderr, "maximal number of ndpoly monomials (1) exceeded!\n");
	  exit(-1);
	}
	monomial_list[i].hbfunction = &ndpoly_heatbath;
	monomial_list[i].accfunction = &ndpoly_acc;
	monomial_list[i].derivativefunction = &ndpoly_derivative;
	no_ndpoly_monomials++;
	monomial_list[i].pf2 = __pf+no*V;
	no++;
      }
      else if(monomial_list[i].type == NDDETRATIO) {
	monomial_list[i].hbfunction = &dummy_heatbath;
	monomial_list[i].accfunction = &nddetratio_acc;
	monomial_list[i].derivativefunction = &dummy_derivative;
	monomial_list[i].pf2 = __pf+no*V;
	monomial_list[i].timescale = -5;
	no++;
      }
    }
    else {
      monomial_list[i].pf = NULL;
      if(no_gauge_monomials > 0) {
	fprintf(stderr, "maximal number of gauge monomials exceeded!\n");
	exit(-1);
      }      
      else if(monomial_list[i].type == GAUGE) {
	monomial_list[i].hbfunction = &gauge_heatbath;
	monomial_list[i].accfunction = &gauge_acc;
	monomial_list[i].derivativefunction = &gauge_derivative;
	no_gauge_monomials++;
	
	if(!monomial_list[i].use_rectangles) {
	  monomial_list[i].c1 = 0.;
	  monomial_list[i].c0 = 1.;
	}
	g_rgi_C1 = monomial_list[i].c1;
	monomial_list[i].c0 = 1. - 8.*monomial_list[i].c1;
	g_rgi_C0 = monomial_list[i].c0;
      }
      else if(monomial_list[i].type == SFGAUGE) {
	monomial_list[i].hbfunction = &sf_gauge_heatbath;
	monomial_list[i].accfunction = &sf_gauge_acc;
	monomial_list[i].derivativefunction = &sf_gauge_derivative;
	no_gauge_monomials++;	
	
	if(!monomial_list[i].use_rectangles) {
	  monomial_list[i].c1 = 0.;
	  monomial_list[i].c0 = 1.;
	  monomial_list[i].c1ss = 0.;
	  monomial_list[i].c1tss = 0.;
	  monomial_list[i].c1tts = 0.;
	}
	g_rgi_C1 = monomial_list[i].c1;
	monomial_list[i].c0 = 1. - 8.*monomial_list[i].c1;
	g_rgi_C0 = monomial_list[i].c0;
	g_C1ss = monomial_list[i].c1ss;
	g_C1tss = monomial_list[i].c1tss;
	g_C1tts = monomial_list[i].c1tts;
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

void dummy_derivative(const int id) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_derivative was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return;
}

void dummy_heatbath(const int id) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_heatbath was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return;
}

double dummy_acc(const int id) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_acc was called. Was that really intended?\n");
    fprintf(stderr, "callers monomial ID was %d\n", id);
  }
  return(0.);
}
