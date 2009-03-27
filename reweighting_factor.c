/***********************************************************************
 * $Id$
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
#include "monomial.h"
#include "reweighting_factor.h"

void reweighting_factor(const int N, const int nstore) {
  int i, j, n = VOLUME;
  double sq_norm, sum=0., x;
  monomial * mnl;
  FILE * ofs;

  for(i = 0; i < N; i++) {
    sq_norm = 0.;
    for(j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	if(mnl->even_odd_flag) {
	  n = VOLUME/2;
	}
	random_spinor_field(mnl->pf, n, mnl->rngrepro);
	mnl->energy0 = square_norm(mnl->pf, n, 1);
      }
    }

    for(j = 0; j < no_monomials; j++) {
      mnl = &monomial_list[j];
      if(mnl->type != GAUGE) {
	sq_norm += -mnl->accfunction(j);
      }
    }

    x = exp(sq_norm);
    sum += x;
    if(g_proc_id == 0 && g_debug_level > 0) {
      printf("rew: sq_norm = %e, W = %e\n", sq_norm, x);
    }
  }
  
  if(g_proc_id == 0) {
    ofs = fopen("reweighing_factor.data", "a");
    fprintf(ofs, "%d %e\n", nstore, sum/N);
    fclose(ofs);
  }
}

