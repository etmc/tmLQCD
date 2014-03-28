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
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "reweighting_factor.h"

void reweighting_factor(const int N, const int nstore) {
  int i, j, n = VOLUME;
  double sq_norm, x, y;
  double * sum, * sum_sq;
  monomial * mnl;
  FILE * ofs;
  hamiltonian_field_t hf;

  hf.gaugefield = g_gauge_field;
  hf.momenta = NULL;
  hf.derivative = NULL;
  hf.update_gauge_copy = g_update_gauge_copy;

  sum = (double*)calloc(no_monomials, sizeof(double));
  sum_sq = (double*)calloc(no_monomials, sizeof(double));

  for(i = 0; i < N; i++) {
    sq_norm = 0.;
    for(j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	if(mnl->even_odd_flag) {
	  random_spinor_field_eo(mnl->pf, mnl->rngrepro, RN_GAUSS);
	}
	else random_spinor_field_lexic(mnl->pf, mnl->rngrepro, RN_GAUSS);
	mnl->energy0 = square_norm(mnl->pf, n, 1);
	if(mnl->type == NDDETRATIO) {
	  if(mnl->even_odd_flag) {
	    random_spinor_field_eo(mnl->pf2, mnl->rngrepro, RN_GAUSS);
	  }
	  else random_spinor_field_lexic(mnl->pf, mnl->rngrepro, RN_GAUSS);
	  mnl->energy0 += square_norm(mnl->pf2, n, 1);
	}
      }
    }

    for(j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	y = mnl->accfunction(j, &hf);
	sq_norm -= y;
	x = exp(sq_norm);
	sum[j] += x;
	sum_sq[j] += x*x;
	if(g_proc_id == 0 && g_debug_level > 0) {
	  printf("monomial[%d] %s, w_%d=%e W=%e\n", j, mnl->name, j, y, x);
	}
      }
    }
  }
  
  if(g_proc_id == 0) {
    ofs = fopen("reweighting_factor.data", "a");
    fprintf(ofs, "%d ", nstore);
    for(j = 0; j < no_monomials; j++) {
      fprintf(ofs, "%e %e ", sum[j]/N, sqrt((-sum[j]*sum[j]/N/N + sum_sq[j]/N)/(N-1)/N));
    }
    fprintf(ofs, "\n");
    fclose(ofs);
  }
}

