/***********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "fatal_error.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"
#include "reweighting_factor.h"

void reweighting_factor(const int N, const int nstore) {
  int n = VOLUME;
  monomial * mnl;
  FILE * ofs;
  hamiltonian_field_t hf;

  hf.gaugefield = g_gauge_field;
  hf.momenta = NULL;
  hf.derivative = NULL;
  hf.update_gauge_copy = g_update_gauge_copy;

  double * data = (double*)calloc(no_monomials*N, sizeof(double));
  double * trlog = (double*)calloc(no_monomials, sizeof(double));

  // we compute the trlog part first, because they are independent of 
  // stochastic noise. This is only needed for even/odd monomials
  for(int j = 0; j < no_monomials; j++) {
    mnl = &monomial_list[j];
    if(mnl->even_odd_flag) {
      init_sw_fields();

      if(mnl->type != NDCLOVERRATCOR && (mnl->kappa != mnl->kappa2
                                       || (mnl->type == NDDETRATIO 
                                           && (mnl->mubar != mnl->mubar2 || mnl->epsbar != mnl->epsbar2))
                                       || (mnl->type != NDDETRATIO
                                           && (mnl->mu != mnl->mu2)))) {
        double c_sw = mnl->c_sw;
        if(c_sw < 0.) c_sw = 0.;
        
        sw_term( (const su3**) hf.gaugefield, mnl->kappa, c_sw); 
        if(mnl->type != NDDETRATIO) {
          trlog[j] = -sw_trace(0, mnl->mu);
        }
        else {
          trlog[j] = -sw_trace_nd(0, mnl->mubar, mnl->epsbar);
        }
        
        sw_term( (const su3**) hf.gaugefield, mnl->kappa2, c_sw);
        if(mnl->type != NDDETRATIO) {
          trlog[j] -= -sw_trace(0, mnl->mu2);
        }
        else {
          trlog[j] -= -sw_trace_nd(0, mnl->mubar2, mnl->epsbar2);
        }
      } else
        trlog[j] = 0.;
    }
    else {
      trlog[j] = 0.;
    }
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("# monomial[%d] %s, trlog = %e\n", j, mnl->name, trlog[j]);
    }
  }

  for(int i = 0; i < N; i++) {
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("# computing reweighting factors for sample %d\n", i);
    }
    for(int j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	if(mnl->even_odd_flag) {
	  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
          mnl->energy0 = square_norm(mnl->pf, n/2, 1);
	}
	else {
          random_spinor_field_lexic(mnl->pf, mnl->rngrepro, RN_GAUSS);
          mnl->energy0 = square_norm(mnl->pf, n, 1);
        }
	if(mnl->type == NDDETRATIO || mnl->type == NDCLOVERRATCOR) {
	  if(mnl->even_odd_flag) {
	    random_spinor_field_eo(mnl->pf2, mnl->rngrepro, RN_GAUSS);
            mnl->energy0 += square_norm(mnl->pf2, n/2, 1);
	  }
	  else {
            random_spinor_field_lexic(mnl->pf2, mnl->rngrepro, RN_GAUSS);
            mnl->energy0 += square_norm(mnl->pf2, n, 1);
          }
	}
	if(g_proc_id == 0 && g_debug_level > 1) {
	  printf("# monomial[%d] %s, energy0 = %e\n", j, mnl->name, mnl->energy0);
	}
      }
    }

    for(int j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	double y = mnl->accfunction(j, &hf);
	data[i*no_monomials + j] = y;
	if(g_proc_id == 0 && g_debug_level > 0) {
	  printf("# monomial[%d] %s, stochastic part: w_%d=%e exp(w_%d)=%e\n", j, mnl->name, j, j, y, exp(y));
	}
      }
    }
  }
  
  if(g_proc_id == 0) {
    char filename[50];
    sprintf(filename, "reweighting_factor.data.%.5d", nstore);
    if((ofs = fopen(filename, "w")) == NULL) {
      fatal_error("Could not open file for data output", "reweighting_factor");
    }
    else {
      for(int j = 0; j < no_monomials; j++) {
        mnl = &monomial_list[j];
        for(int i = 0; i < N; i++) {
          fprintf(ofs, "%.2d %.5d %.12f %.12f %.12f %.12f %.10e\n", j, i, mnl->kappa, mnl->kappa2, mnl->mu, mnl->mu2, data[i*no_monomials + j] + trlog[j]);
        }
      }
      fclose(ofs);
    }
  }
  free(data);
  free(trlog);
}

