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
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "monomial/monomial.h"
#include "xchange/xchange.h"
#include "hamiltonian_field.h"
#include "monitor_forces.h"
#include "gettime.h"

void monitor_forces(hamiltonian_field_t * const hf) {

  for(int id = 0; id < no_monomials; id++) {
    if(monomial_list[ id ].derivativefunction != NULL) {
#ifdef OMP
#pragma omp parallel for
#endif
      for(int i = 0; i < (VOLUMEPLUSRAND + g_dbw2rand);i++) { 
	for(int mu=0;mu<4;mu++) { 
	  _zero_su3adj(hf->derivative[i][mu]);
	}
      }
      
      monomial_list[ id ].derivativefunction(id, hf);
      
#ifdef MPI
      xchange_deri(hf->derivative);
#endif
      
      double sum = 0., max = 0., sum2;
#ifdef OMP
#pragma omp parallel private(sum2)
      {
	int thread_num = omp_get_thread_num();
	g_omp_acc_re[thread_num] = 0.;
#pragma omp for reduction(+ : sum) nowait
#endif
	for(int i = 0; i < VOLUME; i++) {
	  for(int mu = 0; mu < 4; mu++) {
	    sum2 = _su3adj_square_norm(hf->derivative[i][mu]); 
	    sum += sum2;
#ifdef OMP
	    if(sum2 > g_omp_acc_re[thread_num]) g_omp_acc_re[thread_num] = sum2;
#else
	    if(sum2 > max) max = sum2;
#endif
	  }
	}
#ifdef OMP
      } /* OMP closing brace */
      max = g_omp_acc_re[0];
      for( int i = 1; i < omp_num_threads; i++) {
	if(g_omp_acc_re[i] > max) max = g_omp_acc_re[i];
      }
#endif
      
      // output for force monitoring
#ifdef MPI
      MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum = sum2;
      MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max = sum2;
#endif
      if(g_proc_id == 0) {
	printf("# squared force for monomial %s on timescale %d: aver: %1.2e max: %1.2e\n", 
	       monomial_list[ id ].name,
	       monomial_list[ id ].timescale,
	       sum/((double)(VOLUME*g_nproc))/4., max);
	fflush(stdout);
      }
    }
  } 
  return;
}
