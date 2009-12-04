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
#include "observables.h"
#include "hybrid_update.h"
#include "ranlxd.h"
#include "read_input.h"
#include "linsolve.h"
#include "expo.h"
#include "xchange.h"
#include "measure_rectangles.h"
#include "init_gauge_tmp.h"
#include "solver/chrono_guess.h"
#include "solver/bicgstab_complex.h"
#include "update_backward_gauge.h"
#include "update_tm.h"
#include "stout_smear.h"
#include "solver/solver.h"
#include "init_stout_smear_vars.h"
#include "monomial.h"
#include "integrator.h"
/* for the SF: */
#include "sf_calc_action.h"


extern su3 ** g_gauge_field_saved;
void stout_smear();
void unstout();

int update_tm(double *plaquette_energy, double *rectangle_energy, 
	      char * filename, const int return_check, const int acctest) {

  su3 *v, *w;
  static int ini_g_tmp = 0;
  int ix, mu, accept, i=0, j=0;

  double yy[1];
  double dh, expmdh, ret_dh=0., ret_gauge_diff=0., tmp;
  double atime=0., etime=0.;
  double ks,kc,ds,tr,ts,tt;

  /* Energy corresponding to the Gauge part */
  double new_plaquette_energy=0., new_rectangle_energy = 0.;

  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0., ret_enep = 0.;

  /* Energy corresponding to the pseudo fermion part(s) */
  FILE * datafile=NULL, * ret_check_file=NULL;

  /* SF variables: corresponding to the coupling constant */
  double normalisation=0., partial_effective_action=0., coupling=0.; 
  paramsXlfInfo *xlfInfo;

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
      v=&g_gauge_field[ix][mu];
      w=&gauge_tmp[ix][mu];
      _su3_assign(*w,*v);
    }
  }

  /* smear the gauge field */
  if(use_stout_flag == 1) {
    if (bc_flag == 0) {
      stout_smear();
    }
  }
  /* heatbath for all monomials */
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      monomial_list[ Integrator.mnls_per_ts[i][j] ].hbfunction(Integrator.mnls_per_ts[i][j]);
    }
  }

  /* keep on going with the unsmeared gauge field */
  if(use_stout_flag == 1) {
    if (bc_flag == 0) {
      unstout();
    }
  }


  /* initialize the momenta  */
  enep = ini_momenta(reproduce_randomnumber_flag);

  g_sloppy_precision = 1;

  /*run the trajectory*/
  Integrator.integrate[Integrator.no_timescales-1](Integrator.tau, 
						   Integrator.no_timescales-1, 1);

  g_sloppy_precision = 0;
  /*   smear the gauge field */
  if(use_stout_flag == 1) {
    if (bc_flag == 0) {
      stout_smear();
    }
  }

  /* compute the final energy contributions for all monomials */
  dh = 0.;
  for(i = 0; i < Integrator.no_timescales; i++) {
    for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
      dh += monomial_list[ Integrator.mnls_per_ts[i][j] ].accfunction(Integrator.mnls_per_ts[i][j]);
    }
  }
  
  /*   keep on going with the unsmeared gauge field */
  if(use_stout_flag == 1) {
    if (bc_flag == 0) {
      unstout();
    }
  }

  enepx = moment_energy();

  if (bc_flag == 0) { /* if PBC */
    new_plaquette_energy = measure_gauge_action();
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      new_rectangle_energy = measure_rectangles();
    }
  }
  else if (bc_flag == 1) { /* if SF bc */
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      new_plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_iwasaki(g_Tbsf, g_Cs, g_Ct, g_rgi_C0);
      new_rectangle_energy = (1./(2.*3.))*measure_rectangle_sf_iwasaki(g_Tbsf, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
    }
    else {
      new_plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_weights_improvement(g_Tbsf, g_Cs, g_Ct);
    }
  }
  /* Compute the energy difference */
  dh = dh + (enepx - enep);
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

  if(expmdh > yy[0]) {
    accept = 1;
  }
  else {
    accept = 0;
  }
  if(!acctest) {
    accept = 1;
  }
  /* Here a reversibility test is performed */
  /* The trajectory is integrated back      */
  if(return_check == 1) {
    if(accept == 1) {
      xlfInfo = construct_paramsXlfInfo((*plaquette_energy)/(6.*VOLUME*g_nproc), -1);
      write_gauge_field( "conf.save", 64, xlfInfo);
      free(xlfInfo);
    }
    g_sloppy_precision = 1;
    /* run the trajectory back */
    Integrator.integrate[Integrator.no_timescales-1](-Integrator.tau, 
						     Integrator.no_timescales-1, 1);

    g_sloppy_precision = 0;

    /*   compute the energy contributions from the pseudo-fermions  */
    if(use_stout_flag == 1) {
      if (bc_flag == 0) {
	stout_smear();
      }
    }
    
    ret_dh = 0.;
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
	ret_dh += monomial_list[ Integrator.mnls_per_ts[i][j] ].accfunction(Integrator.mnls_per_ts[i][j]);
      }
    }

    /*   keep on going with the unsmeared gauge field */
    if(use_stout_flag == 1) {
      if (bc_flag == 0) {
	unstout();
      }
    }

    ret_enep = moment_energy();

    /* Compute the energy difference */
    ret_dh += (ret_enep - enep );

    /* Compute Differences in the fields */
    ks = 0.;
    kc = 0.;

    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
        tmp = 0.;
        v=&g_gauge_field[ix][mu];
        w=&gauge_tmp[ix][mu];
        ds = ((*v).c00.re-(*w).c00.re)*((*v).c00.re-(*w).c00.re)
          + ((*v).c00.im-(*w).c00.im)*((*v).c00.im-(*w).c00.im)
          + ((*v).c01.re-(*w).c01.re)*((*v).c01.re-(*w).c01.re)
          + ((*v).c01.im-(*w).c01.im)*((*v).c01.im-(*w).c01.im)
          + ((*v).c02.re-(*w).c02.re)*((*v).c02.re-(*w).c02.re)
          + ((*v).c02.im-(*w).c02.im)*((*v).c02.im-(*w).c02.im)
          + ((*v).c10.re-(*w).c10.re)*((*v).c10.re-(*w).c10.re)
          + ((*v).c10.im-(*w).c10.im)*((*v).c10.im-(*w).c10.im)
          + ((*v).c11.re-(*w).c11.re)*((*v).c11.re-(*w).c11.re)
          + ((*v).c11.im-(*w).c11.im)*((*v).c11.im-(*w).c11.im)
          + ((*v).c12.re-(*w).c12.re)*((*v).c12.re-(*w).c12.re)
          + ((*v).c12.im-(*w).c12.im)*((*v).c12.im-(*w).c12.im)
          + ((*v).c20.re-(*w).c20.re)*((*v).c20.re-(*w).c20.re)
          + ((*v).c20.im-(*w).c20.im)*((*v).c20.im-(*w).c20.im)
          + ((*v).c21.re-(*w).c21.re)*((*v).c21.re-(*w).c21.re)
          + ((*v).c21.im-(*w).c21.im)*((*v).c21.im-(*w).c21.im)
          + ((*v).c22.re-(*w).c22.re)*((*v).c22.re-(*w).c22.re)
          + ((*v).c22.im-(*w).c22.im)*((*v).c22.im-(*w).c22.im);
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

    if(accept == 1) {
      read_gauge_field( "conf.save", NULL, NULL, NULL);
    }
  } /* end of reversibility check */

  if(accept == 1) {
    /* accept */
    (*plaquette_energy)=new_plaquette_energy;
    (*rectangle_energy)=new_rectangle_energy;
    /* put the links back to SU(3) group */
    if (bc_flag == 0) { /* if PBC */
      for(ix=0;ix<VOLUME;ix++) { 
	for(mu=0;mu<4;mu++) { 
	  v=&g_gauge_field[ix][mu];
	  *v=restoresu3(*v); 
	}
      }      
    }    
    else if (bc_flag == 1) { /* if SF bc */
      for(ix=0;ix<VOLUME;ix++) { 
	for(mu=0;mu<4;mu++) {
	  if (g_t[ix] == 0 && mu != 0) {
	    v=&g_gauge_field[ix][mu];
	    /* here we do not need to 'restoresu3' because these links are constant ==> we do not want to change them */
	  }
	  else if (g_t[ix]  == g_Tbsf) {
	    v=&g_gauge_field[ix][mu];
	    /* here we do not need to 'restoresu3' because of two reasons: these links are
	       1) either zero ==> they keep updating to zero value all time
	       2) or constant ==> we do not want to change them */
	  }
	  else {
	    v=&g_gauge_field[ix][mu];
	    /* the next line: keeps unitary the gauge field which has been updated */
	    *v=restoresu3(*v);
	  }
	}
      }
    }
  }
  else {
    /* reject: copy gauge_tmp to g_gauge_field */
    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
        v=&g_gauge_field[ix][mu];
        w=&gauge_tmp[ix][mu];
        _su3_assign(*v,*w);
      }
    }
  }
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
#ifdef MPI
  xchange_gauge();
  etime = MPI_Wtime();
#else
  etime = (double)clock()/((double)(CLOCKS_PER_SEC));
#endif

  /* printining data in the .data file */
  if(g_proc_id==0) {
    datafile = fopen(filename, "a");
    if (bc_flag == 0) { /* if PBC */
      fprintf(datafile, "%14.12f %14.12f %e ",
	      (*plaquette_energy)/(6.*VOLUME*g_nproc), dh, expmdh);
    }
    else if (bc_flag == 1) { /* if SF */
      if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) { /* SF with rectangle working */
	
	normalisation = partial_lattice_lo_effective_iwasaki_action_sf_k(g_Tbsf, g_beta, g_rgi_C0, g_rgi_C1, g_eta);
	partial_effective_action = partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
	coupling = normalisation/partial_effective_action;
	
        fprintf(datafile,"%14.12f %14.12f %14.12f %14.12f %14.12f %e ", coupling, normalisation , partial_effective_action, (*plaquette_energy)/(6.*VOLUME*g_nproc), dh, expmdh);
      }
      else { /* SF with only Wilson */

	normalisation = partial_lattice_lo_effective_plaquette_action_sf_k(g_Tbsf, g_beta, g_Ct, g_eta);
	partial_effective_action = partial_wilson_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct);
	coupling = normalisation/partial_effective_action;
	
	fprintf(datafile,"%14.12f %14.12f %14.12f %14.12f %14.12f %e ",
		coupling, normalisation, partial_effective_action, (*plaquette_energy)/(6.*VOLUME*g_nproc), dh, expmdh);
      }
    }
    for(i = 0; i < Integrator.no_timescales; i++) {
      for(j = 0; j < Integrator.no_mnls_per_ts[i]; j++) {
	if(monomial_list[ Integrator.mnls_per_ts[i][j] ].type != GAUGE
	   && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != SFGAUGE 
	   && monomial_list[ Integrator.mnls_per_ts[i][j] ].type != NDPOLY) {
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

void stout_smear() {
  int ix, mu;
    for(ix = 0; ix < VOLUME; ix++) {
      for(mu = 0; mu < 4; mu++) {
	_su3_assign(g_gauge_field_saved[ix][mu], g_gauge_field[ix][mu]);
      }
      stout_smear_gauge_field(stout_rho , stout_no_iter);
    }

  return;
}
void unstout() {
  int ix, mu;
  for(ix = 0; ix < VOLUME; ix++) {
    for(mu = 0; mu < 4; mu++) {
      _su3_assign(g_gauge_field[ix][mu], g_gauge_field_saved[ix][mu]);
    }
  }
  g_update_gauge_copy = 1;
  return;
}


static char const rcsid[] = "$Id$";
