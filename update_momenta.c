/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasebusch
 *               2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#include "update_momenta.h"
#include "clover_leaf.h"
#include "read_input.h"

/* Updates the momenta: equation 16 of Gottlieb */
void update_momenta(int * mnllist, double step, const int no) {
  int i,mu, k;
  double tmp;
  su3adj *xm,*deriv;
  double sum=0., max=0.;
  double sum2=0.;
  double atime=0., etime=0.;

#ifdef MPI
    atime = MPI_Wtime();
#else
    atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

  for(i=0;i<(VOLUME);i++) { 
    for(mu=0;mu<4;mu++) { 
      _zero_su3adj(df0[i][mu]);
    }
  }

  for(k = 0; k < no; k++) {
    /* these are needed for the clover term */
    if(monomial_list[ mnllist[k] ].c_sw > 0) {
      for(i = 0; i < VOLUME; i++) { 
	for(mu = 0; mu < 4; mu++) { 
	  _su3_zero(swm[i][mu]); 
	  _su3_zero(swp[i][mu]); 
	}
      }
    }

    sum = 0.;
    max = 0.;
    for(i = (VOLUME); i < (VOLUME+RAND); i++) { 
      for(mu = 0; mu < 4; mu++) { 
	_zero_su3adj(df0[i][mu]);
      }
    }

    monomial_list[ mnllist[k] ].derivativefunction(mnllist[k]);
#ifdef MPI
    xchange_deri();
#endif
    for(i = 0; i < VOLUME; i++) {
      for(mu = 0; mu < 4; mu++) {
	xm=&moment[i][mu];
	deriv=&df0[i][mu];
	/* force monitoring */
	if(g_debug_level > 0) {
	  sum2 = _su3adj_square_norm(*deriv); 
	  sum+= sum2;
	  if(fabs(sum2) > max) max = sum2;
	}
	tmp = step*monomial_list[ mnllist[k] ].forcefactor;
	/* the minus comes from an extra minus in trace_lambda */
	_minus_const_times_mom(*xm,tmp,*deriv); 
	/* set to zero immediately */
	_zero_su3adj(df0[i][mu]);
      }
    }
#ifdef MPI
    etime = MPI_Wtime();
#else
    etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
    if(g_debug_level > 0) {
#ifdef MPI
      MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum = sum2;
      MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max = sum2;
#endif
      if(g_proc_id == 0) {
	printf("force %s ts %d\taver: %1.4e\tmax: %1.4e\tdt*aver: %1.4e\tdt*max: %1.4e\tt/s %1.4e\n", 
	       monomial_list[ mnllist[k] ].name, monomial_list[ mnllist[k] ].timescale, 
	       fabs(sum/((double)(VOLUME*g_nproc))/4.*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(max*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(step*sum/((double)(VOLUME*g_nproc))/4.*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(step*max*monomial_list[ mnllist[k] ].forcefactor), 
	       etime-atime);
	fflush(stdout);
      }
    }
  }
  return;
}

