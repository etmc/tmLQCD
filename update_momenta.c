/* $Id$ */

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
#include "update_momenta.h"


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
	/* This 2* is coming from what?             */
	/* From a missing factor 2 in trace_lambda? */
	tmp = 2.*step*monomial_list[ mnllist[k] ].forcefactor;
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
	       fabs(2.*sum/((double)(VOLUME*g_nproc))/4.*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(2.*max*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(2.*step*sum/((double)(VOLUME*g_nproc))/4.*monomial_list[ mnllist[k] ].forcefactor),
	       fabs(2.*step*max*monomial_list[ mnllist[k] ].forcefactor), 
	       etime-atime);
	fflush(stdout);
      }
    }
  }
  return;
}

