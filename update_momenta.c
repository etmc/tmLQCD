/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *               2002,2003,2004,2005,2006,2007,2008,2012 Carsten Urbach
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
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "monomial.h"
#include "xchange_deri.h"
#include "clover_leaf.h"
#include "read_input.h"
#include "hamiltonian_field.h"
#include "update_momenta.h"

/* Updates the momenta: equation 16 of Gottlieb */
void update_momenta(int * mnllist, double step, const int no, 
		    hamiltonian_field_t * const hf) {
  int i,mu, k;
  su3adj *xm,*deriv;
  double atime=0., etime=0.;

  for(i = 0; i < (VOLUMEPLUSRAND);i++) { 
    for(mu=0;mu<4;mu++) { 
      _zero_su3adj(hf->derivative[i][mu]);
    }
  }

  for(k = 0; k < no; k++) {
    if(monomial_list[ mnllist[k] ].derivativefunction != NULL) {
#ifdef MPI
      atime = MPI_Wtime();
#else
      atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

      /* these are needed for the clover term */
      if(monomial_list[ mnllist[k] ].type == 9 || monomial_list[ mnllist[k] ].type == 10) {
	for(i = 0; i < VOLUME; i++) { 
	  for(mu = 0; mu < 4; mu++) { 
	    _su3_zero(swm[i][mu]);
	    _su3_zero(swp[i][mu]);
	  }
	}
      }
      
      monomial_list[ mnllist[k] ].derivativefunction(mnllist[k], hf);
#ifdef MPI
      etime = MPI_Wtime();
#else
      etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
    }
  }
#ifdef MPI
  xchange_deri(hf->derivative);
#endif
  for(i = 0; i < VOLUME; i++) {
    for(mu = 0; mu < 4; mu++) {
      xm=&hf->momenta[i][mu];
      deriv=&hf->derivative[i][mu];
      /* the minus comes from an extra minus in trace_lambda */
      _su3adj_minus_const_times_su3adj(*xm, step, *deriv); 
    }
  }
  return;
}

