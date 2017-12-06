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
 *
 * This routine contains the update part for
 * the HMC with up to three pseudo fermion fields
 * for twisted mass QCD
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *
 * Modified by Jenifer Gonzalez Lopez for the Schroedinger Functional
 *
 ***********************************************************************/

#include <lime.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef TM_USE_MPI
# include <mpi.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
#endif
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "operator/tm_operators.h"
#include "linalg_eo.h"
#include "io/gauge.h"
#include "io/params.h"
#include "measure_gauge_action.h"
#include "ranlxd.h"
#include "read_input.h"
#include "expo.h"
#include "xchange/xchange.h"
#include "measure_rectangles.h"
#include "init/init_gauge_tmp.h"
#include "init/init_gauge_field.h"
#include "monomial/monomial.h"
#include "integrator.h"
#include "hamiltonian_field.h"
#include "update_smd.h"
#include "gettime.h"
#ifdef DDalphaAMG
#include "DDalphaAMG_interface.h"
#endif

extern su3 ** g_gauge_field_saved;

void update_smd(double *plaquette_energy, double *rectangle_energy, 
                char * filename, const int return_check, const int acctest, 
                const int traj_counter) {
  
  su3 *v, *w;
  int i=0, j=0;

  double dh, expmdh;
  double atime=0., etime=0.;

  /* Energy corresponding to the Gauge part */
  double new_plaquette_energy=0., new_rectangle_energy = 0.;

  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0.;

  /* Energy corresponding to the pseudo fermion part(s) */
  FILE * datafile=NULL;
  hamiltonian_field_t hf;

  hf.gaugefield = g_gauge_field;
  hf.momenta = moment;
  hf.derivative = df0;
  hf.update_gauge_copy = g_update_gauge_copy;
  hf.traj_counter = traj_counter;
  integrator_set_fields(&hf);

  atime = gettime();

  /*
   *  here the momentum and spinor fields are initialized 
   *  and their respective actions are calculated
   */

  /* 
   *  copy the gauge field to gauge_tmp 
   */
#ifdef TM_USE_OMP
#pragma omp parallel for private(w,v)
#endif
  for(int ix=0;ix<VOLUME;ix++) { 
    for(int mu=0;mu<4;mu++) {
      v=&hf.gaugefield[ix][mu];
      w=&gauge_tmp[ix][mu];
      _su3_assign(*w,*v);
    }
  }

#ifdef DDalphaAMG
  MG_reset();
#endif

  /* heatbath for all monomials */
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      monomial_list[ Integrator.mnls_per_ts[i][j] ].hbfunction(Integrator.mnls_per_ts[i][j], &hf);
    }
  }

  if(Integrator.monitor_forces) monitor_forces(&hf);
  /* initialize the momenta  */
  double epsilon = Integrator.tau/(double)Integrator.n_int[Integrator.no_timescales-1];
  enep = random_su3adj_field_smd(reproduce_randomnumber_flag, hf.momenta, epsilon, smd_gamma);
  
  g_sloppy_precision = 1;

  // run the trajectory
  if(Integrator.n_int[Integrator.no_timescales-1] > 0) {
    Integrator.integrate[Integrator.no_timescales-1](Integrator.tau, 
                                                     Integrator.no_timescales-1, 1,
                                                     Integrator.tau);
  }

  g_sloppy_precision = 0;

  // compute the final energy contributions for all monomials
  // this will not be needed once we don't do the accept reject anymore.
  dh = 0.;
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      dh += monomial_list[ Integrator.mnls_per_ts[i][j] ].accfunction(Integrator.mnls_per_ts[i][j], &hf);
    }
  }
  // not needed anymore once we don't do the accept reject anymore.
  enepx = moment_energy(hf.momenta);

  new_plaquette_energy = measure_plaquette( (const su3**) hf.gaugefield);
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    new_rectangle_energy = measure_rectangles( (const su3**) hf.gaugefield);
  }

  if(g_proc_id == 0 && g_debug_level > 3) printf("called moment_energy: dh = %1.10e\n", (enepx - enep));
  /* Compute the energy difference */
  dh = dh + (enepx - enep);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called momenta_acc dH = %e\n", (enepx - enep));
  }

  *plaquette_energy = new_plaquette_energy;
  *rectangle_energy = new_rectangle_energy;
  /* put the links back to SU(3) group */
#ifdef TM_USE_OMP
#pragma omp parallel for private(v)
#endif
  for(int ix=0;ix<VOLUME;ix++) { 
    for(int mu=0;mu<4;mu++) { 
      v=&hf.gaugefield[ix][mu];
      restoresu3_in_place(v); 
    }
  }

  hf.update_gauge_copy = 1;
  g_update_gauge_copy = 1;
  g_update_gauge_copy_32 = 1;  
#ifdef TM_USE_MPI
  xchange_gauge(hf.gaugefield);
#endif
  
  /*Convert to a 32 bit gauge field, after xchange*/
  convert_32_gauge_field(g_gauge_field_32, hf.gaugefield, VOLUMEPLUSRAND + g_dbw2rand); 
  
  etime=gettime();

  /* printing data in the .data file */
  if(g_proc_id==0) {
    datafile = fopen(filename, "a");
    if (!bc_flag) { /* if Periodic Boundary Conditions */
      fprintf(datafile, "%.8d %14.12f %14.12f %e ", traj_counter,
              (*plaquette_energy)/(6.*VOLUME*g_nproc), dh, expmdh);
    }
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
        if(monomial_list[ Integrator.mnls_per_ts[i][j] ].type != GAUGE
            && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != SFGAUGE 
            && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != NDPOLY
            && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != NDCLOVER
            && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != CLOVERNDTRLOG
            && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != CLOVERTRLOG ) {
          fprintf(datafile,"%d %d ",  monomial_list[ Integrator.mnls_per_ts[i][j] ].iter0, 
                  monomial_list[ Integrator.mnls_per_ts[i][j] ].iter1);
        }
      }
    }
    fprintf(datafile, "%e", etime-atime);
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0) {
      fprintf(datafile, " %e", (*rectangle_energy)/(12*VOLUME*g_nproc));
    }
    fprintf(datafile, "\n");
    fflush(datafile);
    fclose(datafile);
  }
  return;
}

