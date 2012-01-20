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
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "tm_operators.h"
#include "linalg_eo.h"
#include "io/gauge.h"
#include "io/params.h"
#include "measure_gauge_action.h"
#include "hybrid_update.h"
#include "ranlxd.h"
#include "read_input.h"
#include "expo.h"
#include "xchange.h"
#include "measure_rectangles.h"
#include "init_gauge_tmp.h"
#include "monomial.h"
#include "integrator.h"
#include "hamiltonian_field.h"
#include "update_tm.h"

extern su3 ** g_gauge_field_saved;

int update_tm(double *plaquette_energy, double *rectangle_energy, 
              char * filename, const int return_check, const int acctest) {

  su3 *v, *w;
  static int ini_g_tmp = 0;
  int ix, mu, accept, i=0, j=0, iostatus=0;

  double yy[1];
  double dh, expmdh, ret_dh=0., ret_gauge_diff=0., tmp;
  double atime=0., etime=0.;
  double ks,kc,ds,tr,ts,tt;

  char tmp_filename[50];

  /* Energy corresponding to the Gauge part */
  double new_plaquette_energy=0., new_rectangle_energy = 0.;

  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0., ret_enep = 0.;

  /* Energy corresponding to the pseudo fermion part(s) */
  FILE * datafile=NULL, * ret_check_file=NULL;
  hamiltonian_field_t hf;
  paramsXlfInfo *xlfInfo;

  hf.gaugefield = g_gauge_field;
  hf.momenta = moment;
  hf.derivative = df0;
  hf.update_gauge_copy = g_update_gauge_copy;
  hf.update_gauge_energy = g_update_gauge_energy;
  hf.update_rectangle_energy = g_update_rectangle_energy;
  integrator_set_fields(&hf);

  strcpy(tmp_filename, ".conf.tmp");
  if(ini_g_tmp == 0) {
    ini_g_tmp = init_gauge_tmp(VOLUME);
    if(ini_g_tmp != 0) {
      exit(-1);
    }
    ini_g_tmp = 1;
  }
#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/((double)(CLOCKS_PER_SEC));
#endif

  /*
   *  here the momentum and spinor fields are initialized 
   *  and their respective actions are calculated
   */

  /* 
   *  copy the gauge field to gauge_tmp 
   */
  for(ix=0;ix<VOLUME;ix++) { 
    for(mu=0;mu<4;mu++) {
      v=&hf.gaugefield[ix][mu];
      w=&gauge_tmp[ix][mu];
      _su3_assign(*w,*v);
    }
  }

  /* heatbath for all monomials */
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      monomial_list[ Integrator.mnls_per_ts[i][j] ].hbfunction(Integrator.mnls_per_ts[i][j], &hf);
    }
  }

  /* initialize the momenta  */
  enep = init_momenta(reproduce_randomnumber_flag, hf.momenta);

  g_sloppy_precision = 1;

  /* run the trajectory */
  Integrator.integrate[Integrator.no_timescales-1](Integrator.tau, 
                       Integrator.no_timescales-1, 1);

  g_sloppy_precision = 0;

  /* compute the final energy contributions for all monomials */
  dh = 0.;
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      dh += monomial_list[ Integrator.mnls_per_ts[i][j] ].accfunction(Integrator.mnls_per_ts[i][j], &hf);
    }
  }

  enepx = moment_energy(hf.momenta);

  if (!bc_flag) { /* if PBC */
    new_plaquette_energy = measure_gauge_action(hf.gaugefield);
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      new_rectangle_energy = measure_rectangles(hf.gaugefield);
    }
  }
  /* Compute the energy difference */
  dh = dh + (enepx - enep);
  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("called momenta_acc dH = %e\n", (enepx - enep));
  }
  expmdh = exp(-dh);
  /* the random number is only taken at node zero and then distributed to 
     the other sites */
  if(g_proc_id==0) {
    ranlxd(yy,1);
#ifdef MPI
    for(i = 1; i < g_nproc; i++) {
      MPI_Send(&yy[0], 1, MPI_DOUBLE, i, 31, MPI_COMM_WORLD);
    }
#endif
  }
#ifdef MPI
  else{
    MPI_Recv(&yy[0], 1, MPI_DOUBLE, 0, 31, MPI_COMM_WORLD, &status);
  }
#endif

  accept = (!acctest | (expmdh > yy[0]));
  if(g_proc_id == 0) {
    fprintf(stdout, "# Trajectory is %saccepted.\n", (accept ? "" : "not "));
  }
  /* Here a reversibility test is performed */
  /* The trajectory is integrated back      */
  if(return_check) {
    if(g_proc_id == 0) {
      fprintf(stdout, "# Performing reversibility check.\n");
    }
    if(accept) {
      /* save gauge file to disk before performing reversibility check */
      xlfInfo = construct_paramsXlfInfo((*plaquette_energy)/(6.*VOLUME*g_nproc), -1);
      // Should write this to temporary file first, and then check
      if(g_proc_id == 0 && g_debug_level > 0) {
        fprintf(stdout, "# Writing gauge field to file %s.\n", tmp_filename);
      }
      if((iostatus = write_gauge_field( tmp_filename, 64, xlfInfo) != 0 )) {
        /* Writing failed directly */
        fprintf(stderr, "Error %d while writing gauge field to %s\nAborting...\n", iostatus, tmp_filename);
        exit(-2);
      }
      /* There is double writing of the gauge field, also in hmc_tm.c in this case */
      /* No reading back check needed here, as reading back is done further down */
      if(g_proc_id == 0 && g_debug_level > 0) {
        fprintf(stdout, "# Writing done.\n");
      }
      free(xlfInfo);
    }
    g_sloppy_precision = 1;
    /* run the trajectory back */
    Integrator.integrate[Integrator.no_timescales-1](-Integrator.tau, 
                         Integrator.no_timescales-1, 1);

    g_sloppy_precision = 0;

    /*   compute the energy contributions from the pseudo-fermions  */
    ret_dh = 0.;
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
        ret_dh += monomial_list[ Integrator.mnls_per_ts[i][j] ].accfunction(Integrator.mnls_per_ts[i][j], &hf);
      }
    }

    ret_enep = moment_energy(hf.momenta);

    /* Compute the energy difference */
    ret_dh += ret_enep - enep ;

    /* Compute Differences in the fields */
    ks = 0.;
    kc = 0.;

    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
        tmp = 0.;
        v=&hf.gaugefield[ix][mu];
        w=&gauge_tmp[ix][mu];
        ds = (creal(v->c00)-creal(w->c00))*(creal(v->c00)-creal(w->c00))
          + (cimag(v->c00)-cimag(w->c00))*(cimag(v->c00)-cimag(w->c00))
          + (creal(v->c01)-creal(w->c01))*(creal(v->c01)-creal(w->c01))
          + (cimag(v->c01)-cimag(w->c01))*(cimag(v->c01)-cimag(w->c01))
          + (creal(v->c02)-creal(w->c02))*(creal(v->c02)-creal(w->c02))
          + (cimag(v->c02)-cimag(w->c02))*(cimag(v->c02)-cimag(w->c02))
          + (creal(v->c10)-creal(w->c10))*(creal(v->c10)-creal(w->c10))
          + (cimag(v->c10)-cimag(w->c10))*(cimag(v->c10)-cimag(w->c10))
          + (creal(v->c11)-creal(w->c11))*(creal(v->c11)-creal(w->c11))
          + (cimag(v->c11)-cimag(w->c11))*(cimag(v->c11)-cimag(w->c11))
          + (creal(v->c12)-creal(w->c12))*(creal(v->c12)-creal(w->c12))
          + (cimag(v->c12)-cimag(w->c12))*(cimag(v->c12)-cimag(w->c12))
          + (creal(v->c20)-creal(w->c20))*(creal(v->c20)-creal(w->c20))
          + (cimag(v->c20)-cimag(w->c20))*(cimag(v->c20)-cimag(w->c20))
          + (creal(v->c21)-creal(w->c21))*(creal(v->c21)-creal(w->c21))
          + (cimag(v->c21)-cimag(w->c21))*(cimag(v->c21)-cimag(w->c21))
          + (creal(v->c22)-creal(w->c22))*(creal(v->c22)-creal(w->c22))
          + (cimag(v->c22)-cimag(w->c22))*(cimag(v->c22)-cimag(w->c22));
        ds = sqrt(ds);
        tr = ds + kc;
        ts = tr + ks;
        tt = ts-ks;
        ks = ts;
        kc = tr-tt;
      }
    }
    ret_gauge_diff = ks + kc;
#ifdef MPI
    tmp = ret_gauge_diff;
    MPI_Reduce(&tmp, &ret_gauge_diff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    /* compute the total H */
    tmp = enep;
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
        tmp += monomial_list[ Integrator.mnls_per_ts[i][j] ].energy0;
      }
    }
    /* Output */
    if(g_proc_id == 0) {
      ret_check_file = fopen("return_check.data","a");
      fprintf(ret_check_file,"ddh = %1.4e ddU= %1.4e ddh/H = %1.4e\n",
              ret_dh, ret_gauge_diff/4./((double)(VOLUME*g_nproc))/3., ret_dh/tmp);
      fclose(ret_check_file);
    }

    if(accept) {
      /* Read back gauge field */
      if(g_proc_id == 0 && g_debug_level > 0) {
        fprintf(stdout, "# Trying to read gauge field from file %s.\n", tmp_filename);
      }

      if((iostatus = read_gauge_field(tmp_filename) != 0)) {
        fprintf(stderr, "Error %d while reading gauge field from %s\nAborting...\n", iostatus, tmp_filename);
        exit(-2);
      }
      if(g_proc_id == 0 && g_debug_level > 0) {
        fprintf(stdout, "# Reading done.\n");
      }
    }
    if(g_proc_id == 0) {
      fprintf(stdout, "# Reversibility check done.\n");
    }
  } /* end of reversibility check */

  if(accept) {
    *plaquette_energy = new_plaquette_energy;
    *rectangle_energy = new_rectangle_energy;
    /* put the links back to SU(3) group */
    if (!bc_flag) { /* periodic boundary conditions */
      for(ix=0;ix<VOLUME;ix++) { 
        for(mu=0;mu<4;mu++) { 
          v=&hf.gaugefield[ix][mu];
          *v=restoresu3(*v); 
        }
      }
    }
  }
  else { /* reject: copy gauge_tmp to hf.gaugefield */
    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
        v=&hf.gaugefield[ix][mu];
        w=&gauge_tmp[ix][mu];
        _su3_assign(*v,*w);
      }
    }
  }
  hf.update_gauge_copy = 1;
  g_update_gauge_copy = 1;
  hf.update_gauge_energy = 1;
  g_update_gauge_energy = 1;
  hf.update_rectangle_energy = 1;
  g_update_rectangle_energy = 1;
#ifdef MPI
  xchange_gauge(hf.gaugefield);
  etime = MPI_Wtime();
#else
  etime = (double)clock()/((double)(CLOCKS_PER_SEC));
#endif

  /* printing data in the .data file */
  if(g_proc_id==0) {
    datafile = fopen(filename, "a");
    if (!bc_flag) { /* if Periodic Boundary Conditions */
      fprintf(datafile, "%14.12f %14.12f %e ",
              (*plaquette_energy)/(6.*VOLUME*g_nproc), dh, expmdh);
    }
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
        if(monomial_list[ Integrator.mnls_per_ts[i][j] ].type != GAUGE
	   && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != SFGAUGE 
	   && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != NDPOLY
	   && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != CLOVERTRLOG ) {
          fprintf(datafile,"%d %d ",  monomial_list[ Integrator.mnls_per_ts[i][j] ].iter0, 
                  monomial_list[ Integrator.mnls_per_ts[i][j] ].iter1);
        }
      }
    }
    fprintf(datafile, "%d %e", accept, etime-atime);
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0) {
      fprintf(datafile, " %e", (*rectangle_energy)/(12*VOLUME*g_nproc));
    }
    fprintf(datafile, "\n");
    fflush(datafile);
    fclose(datafile);
  }
  return(accept);
}

