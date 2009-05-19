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
 *
 * Hybrid-Monte-Carlo for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 * Modified by Jenifer Gonzalez Lopez for the Schroedinger Functional
 *
 *******************************************************************************/
#define MAIN_PROGRAM

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "start.h"
#include "observables.h"
#include "measure_rectangles.h"
#ifdef MPI
# include "xchange.h"
#endif
#ifdef HAVE_LIBLEMON
# ifndef MPI
#  error "Parallel IO without MPI not supported"
# endif
# include "parallel_io.h"
#endif
#include "io.h"
#include "propagator_io.h"
#include "gauge_io.h"
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "update_tm.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_gauge_tmp.h"
#include "init_dirac_halfspinor.h"
#include "init_stout_smear_vars.h"
#include "init_bispinor_field.h"
#include "init_chi_spinor_field.h"
#include "xchange_halffield.h"
#include "test/check_geometry.h"
#include "boundary.h"
#include "phmc.h"
#include "solver/solver.h"
#include "polyakov_loop.h"
#include "monomial.h"
#include "integrator.h"
#include "sighandler.h"
#include "online_measurement.h"
#include "sf_calc_action.h"
#include "sf_observables.h"

void usage(){
  fprintf(stdout, "HMC for Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   hmc_tm [options]\n");
  fprintf(stdout, "Options: [-f input-filename]  default: hmc.input\n");
  fprintf(stdout, "         [-o output-filename] default: output\n");
  fprintf(stdout, "         [-v] more verbosity\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}

extern int nstore;

int const rlxdsize = 105;

int main(int argc,char *argv[]) {

  FILE *parameterfile=NULL, *countfile=NULL;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char gauge_filename[50];
  char nstore_filename[50];
  char tmp_filename[50];
  char * input_filename = NULL;

  int rlxd_state[rlxdsize];

  int j,ix,mu, trajectory_counter=1;
  struct timeval t1;

  /* Energy corresponding to the Gauge part */
  double eneg = 0., plaquette_energy = 0., rectangle_energy = 0.;
  /* Acceptance rate */
  int Rate=0;
  /* Do we want to perform reversibility checks */
  /* See also return_check_flag in read_input.h */
  int return_check = 0;
  /* For getopt */
  int c;

  /* For the Polyakov loop: */
  int dir = 2;
  complex pl, pl4;
#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL,&catch_ill_inst);
#endif

  strcpy(gauge_filename,"conf.save");
  strcpy(nstore_filename,".nstore_counter");
  strcpy(tmp_filename, ".conf.tmp");

  verbose = 0;
  g_use_clover_flag = 0;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);
#else
  g_proc_id = 0;
#endif


  while ((c = getopt(argc, argv, "h?vf:o:")) != -1) {
    switch (c) {
    case 'f':
      input_filename = calloc(200, sizeof(char));
      strcpy(input_filename,optarg);
      break;
    case 'o':
      filename = calloc(200, sizeof(char));
      strcpy(filename,optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
    case '?':
    default:
      usage();
      break;
    }
  }
  if(input_filename == NULL){
    input_filename = "hmc.input";
  }
  if(filename == NULL){
    filename = "output";
  }

  /* Read the input file */
  read_input(input_filename);

  DUM_DERI = 6;
  DUM_SOLVER = DUM_DERI+8;
  DUM_MATRIX = DUM_SOLVER+6;
  if(g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX+8;
  }
  else {
    NO_OF_SPINORFIELDS = DUM_MATRIX+6;
  }
  DUM_BI_DERI = 6;
  DUM_BI_SOLVER = DUM_BI_DERI+7;
  DUM_BI_MATRIX = DUM_BI_SOLVER+6;
  NO_OF_BISPINORFIELDS = DUM_BI_MATRIX+6;

  mpi_init(argc, argv);

  init_integrator();

  if(nstore == -1) {
    countfile = fopen(nstore_filename, "r");
    if(countfile != NULL) {
      j = fscanf(countfile, "%d %d %s\n", &nstore, &trajectory_counter, gauge_input_filename);
      if(j < 1) nstore = 0;
      if(j < 2) trajectory_counter = 0;
      fclose(countfile);
    }
    else {
      nstore = 0;
      trajectory_counter = 0;
    }
  }

  if(g_rgi_C1 == 0.) {
    g_dbw2rand = 0;
  }
#ifndef MPI
  g_dbw2rand = 0;
#endif

  if(g_proc_id == 0) {
    for(j = 0; j < no_monomials; j++) {
      printf("# monomial id %d type = %d timescale %d\n", j, monomial_list[j].type, monomial_list[j].timescale);
    }
  }

  g_mu = g_mu1;

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(0);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry_indices! Aborting...\n");
    exit(0);
  }
  if(even_odd_flag) {
    j = init_monomials(VOLUMEPLUSRAND/2, even_odd_flag);
  }
  else {
    j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for monomial pseudo fermion  fields! Aborting...\n");
    exit(0);
  }
  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  if(even_odd_flag) {
    j = init_csg_field(VOLUMEPLUSRAND/2);
  }
  else {
    j = init_csg_field(VOLUMEPLUSRAND);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for csg fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }

  if(g_running_phmc) {
    j = init_bispinor_field(VOLUME/2, NO_OF_BISPINORFIELDS);
    if (j!= 0) {
      fprintf(stderr, "Not enough memory for Bispinor fields! Aborting...\n");
      exit(0);
    }
  }

  zero_spinor_field(g_spinor_field[DUM_DERI+4],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+5],VOLUME);
  zero_spinor_field(g_spinor_field[DUM_DERI+6],VOLUME);

  if(use_stout_flag == 1)
    init_stout_smear_vars(VOLUMEPLUSRAND, stout_no_iter);

  /*construct the filenames for the observables and the parameters*/
  strcpy(datafilename,filename);  strcat(datafilename,".data");
  strcpy(parameterfilename,filename);  strcat(parameterfilename,".para");

  if(g_proc_id == 0){
    parameterfile = fopen(parameterfilename, "a");
    write_first_messages(parameterfile, 0);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  check_geometry();

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    init_dirac_halfspinor32();
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif

  /* Initialise random number generator */
  start_ranlux(rlxd_level, random_seed^nstore);

  /* Set up the gauge field */
  /* continue and restart */
  if(startoption==3 || startoption == 2) {
    if(g_proc_id == 0) {
      printf("# Reading Gauge field from file %s in %d Bit\n",
	     gauge_input_filename, gauge_precision_read_flag);
      fflush(stdout);
    }
#ifdef HAVE_LIBLEMON
      read_lemon_gauge_field_parallel(gauge_input_filename, rlxd_state);
#else
    if(gauge_precision_read_flag == 64) {
      read_lime_gauge_field(gauge_input_filename);
    }
    else if(gauge_precision_read_flag == 32){
      read_lime_gauge_field_singleprec(gauge_input_filename);
    }
    #endif

    if (g_proc_id == 0){
      printf("# done!\n"); fflush(stdout);
    }
  }
  else if (startoption == 1) {
    /* hot */
    random_gauge_field(reproduce_randomnumber_flag);
  }
  else if(startoption == 0) {
    /* cold */
    unit_g_gauge_field();
  }

  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge();
#endif

  if(g_running_phmc) init_phmc();


  /*********************************************************/
  /* impose SF bc in case it was chosen in the input file */
  /*******************************************************/

  if (bc_flag == 1) { /* if SF */
    dirichlet_boundary_conditions(g_Tbsf);
    dirichlet_boundary_conditions_spatial_links_to_one(g_Tbsf);
    //sf_boundary_conditions_spatially_constant_abelian_field(g_Tbsf, g_eta);
    fprintf(parameterfile,"# SF put boundary at time slice: g_Tbsf = %d \n",g_Tbsf);

    /* compute the energy of the gauge field for SF */
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      /* NOTE: the factor (1./(2.*3.)) is due to the difference between	our normalisation and Carstens's normalisation
	 when defining the plaquette and rectangle functions */
      plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_iwasaki(g_Tbsf, g_Cs, g_Ct, g_rgi_C0);
      rectangle_energy = (1./(2.*3.))*measure_rectangle_sf_iwasaki(g_Tbsf, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
      eneg = plaquette_energy + rectangle_energy;
      /* print energy for SF */
      if(g_proc_id==0){
	fprintf(parameterfile,"# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	fprintf(parameterfile,"# First rectangle value for SF: %14.12f \n", rectangle_energy/(12.*VOLUME*g_nproc));
	printf("# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	printf("# First rectangle value for SF: %14.12f \n", rectangle_energy/(12.*VOLUME*g_nproc));
      }
    }
    else {
      /* NOTE: the factor (1./(2.*3.)) is due to the difference between	our normalisation and Carstens's normalisation
	 when defining the plaquette and rectangle functions */
      plaquette_energy = (1./(2.*3.))*measure_plaquette_sf_weights_improvement(g_Tbsf, g_Cs, g_Ct);
      eneg = plaquette_energy;
      /* print plaquette energy for SF */
      if(g_proc_id==0){
	fprintf(parameterfile,"# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
	printf("# First plaquette value for SF: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      }
    }
  }
  else if (bc_flag == 0) { /*if PBC */
    plaquette_energy=measure_gauge_action();
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      rectangle_energy = measure_rectangles();
      if(g_proc_id==0){
	fprintf(parameterfile,"# First rectangle value: %14.12f \n",rectangle_energy/(12.*VOLUME*g_nproc));
      }
    }
    eneg = g_rgi_C0 * plaquette_energy + g_rgi_C1 * rectangle_energy;

    /* Measure and print the Polyakov loop: */
    polyakov_loop(&pl, dir);

    if(g_proc_id==0){
      fprintf(parameterfile,"# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      printf("# First plaquette value: %14.12f \n", plaquette_energy/(6.*VOLUME*g_nproc));
      fprintf(parameterfile,"# First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	      dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
    }

    dir=3;
    polyakov_loop(&pl, dir);
    if(g_proc_id==0){
      fprintf(parameterfile,"# First Polyakov loop value in %d-direction |L(%d)|= %14.12f \n",
	      dir, dir, sqrt(pl.re*pl.re+pl.im*pl.im));
      fclose(parameterfile);
    }
  }


  /* set ddummy to zero */
  for(ix = 0; ix < VOLUME+RAND; ix++){
    for(mu=0; mu<4; mu++){
      ddummy[ix][mu].d1=0.;
      ddummy[ix][mu].d2=0.;
      ddummy[ix][mu].d3=0.;
      ddummy[ix][mu].d4=0.;
      ddummy[ix][mu].d5=0.;
      ddummy[ix][mu].d6=0.;
      ddummy[ix][mu].d7=0.;
      ddummy[ix][mu].d8=0.;
    }
  }

  if(g_proc_id == 0) {
    gettimeofday(&t1,NULL);
    countfile = fopen("history_hmc_tm", "a");
    fprintf(countfile, "!!! Timestamp %ld, Nsave = %d, g_mu = %e, g_mu1 = %e, g_mu_2 = %e, g_mu3 = %e, beta = %f, kappa = %f, C1 = %f, ",
	    t1.tv_sec, Nsave, g_mu, g_mu1, g_mu2, g_mu3, g_beta, g_kappa, g_rgi_C1);
    for(j = 0; j < Integrator.no_timescales; j++) {
      fprintf(countfile, "n_int[%d] = %d ", j, Integrator.no_mnls_per_ts[j]);
    }
    fprintf(countfile, "\n");
    fclose(countfile);
  }

  /* Loop for measurements */
  for(j = 0; j < Nmeas; j++) {
    if(g_proc_id == 0) {
      printf("# Starting trajectory no %d\n", trajectory_counter);
    }

    if(return_check_flag == 1 && trajectory_counter%return_check_interval == 0) return_check = 1;
    else return_check = 0;

    Rate += update_tm(&plaquette_energy, &rectangle_energy, datafilename, return_check, Ntherm<trajectory_counter);
    /*     Rate += update_tm(integtyp, &plaquette_energy, &rectangle_energy, datafilename,  */
    /* 		      dtau, Nsteps, nsmall, tau, int_n, return_check, lambda, reproduce_randomnumber_flag); */


    if (bc_flag == 0) { /* if PBC */
      /* Measure the Polyakov loop in direction 2 and 3:*/
      polyakov_loop(&pl, 2);
      polyakov_loop(&pl4, 3);
    }

    /* Save gauge configuration all Nsave times */
    if((Nsave !=0) && (trajectory_counter%Nsave == 0) && (trajectory_counter!=0)) {
      sprintf(gauge_filename,"conf.%.4d", nstore);
      if(g_proc_id == 0) {
        countfile = fopen("history_hmc_tm", "a");
	fprintf(countfile, "%.4d, measurement %d of %d, Nsave = %d, Plaquette = %e, |L(%d)| = %e, |L(%d)| = %e trajectory nr = %d\n",
		nstore, j, Nmeas, Nsave, plaquette_energy/(6.*VOLUME*g_nproc),
		2, sqrt(pl.re*pl.re+pl.im*pl.im),
		dir, sqrt(pl4.re*pl4.re+pl4.im*pl4.im), trajectory_counter);
	fclose(countfile);
      }
      nstore ++;
    }
    else {
      sprintf(gauge_filename,"conf.save");
    }

    if(((Nsave !=0) && (trajectory_counter%Nsave == 0) && (trajectory_counter!=0)) || (write_cp_flag == 1) || (j >= (Nmeas - 1))) {
      /* Write the gauge configuration first to a temporary file */
/*       write_gauge_field_time_p( tmp_filename); */
#ifdef HAVE_LIBLEMON
      write_lemon_gauge_field_parallel( tmp_filename , plaquette_energy/(6.*VOLUME*g_nproc),
                  trajectory_counter, gauge_precision_write_flag);
#else
      write_lime_gauge_field( tmp_filename , plaquette_energy/(6.*VOLUME*g_nproc),
                  trajectory_counter, gauge_precision_write_flag);
      /*  write the status of the random number generator to a file */
      if(g_proc_id==0) {
    rlxd_get(rlxd_state);
    write_rlxd_state(tmp_filename, rlxd_state, rlxdsize);
      }
#endif


      /* Now move it! */
      if(g_proc_id == 0) {
  	rename(tmp_filename, gauge_filename);
        countfile = fopen(nstore_filename, "w");
        fprintf(countfile, "%d %d %s\n", nstore, trajectory_counter+1, gauge_filename);
        fclose(countfile);
      }
    }
    if(online_measurement_flag && (trajectory_counter%online_measurement_freq == 0)) {
      online_measurement(trajectory_counter,
			 ((int)(100000*plaquette_energy/(6.*VOLUME*g_nproc)))%(g_nproc_t*T));
    }

    if((g_rec_ev !=0) && (trajectory_counter%g_rec_ev == 0) && (g_running_phmc)) {
      phmc_compute_ev(trajectory_counter, plaquette_energy);
    }


    if(g_proc_id == 0) {
      verbose = 1;
    }
    ix = reread_input("hmc.reread");
    if(g_proc_id == 0) {
      verbose = 0;
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(ix == 0 && g_proc_id == 0) {
      countfile = fopen("history_hmc_tm", "a");
      fprintf(countfile, "# Changed parameter according to hmc.reread: measurement %d of %d\n", j, Nmeas);
      fclose(countfile);
      printf("# Changed parameter according to hmc.reread (see stdout): measurement %d of %d\n", j, Nmeas);
      remove("hmc.reread");
    }
    trajectory_counter++;
  } /* end of loop over trajectories */

  if(g_proc_id==0) {
    printf("Acceptance Rate was: %e Percent\n", 100.*(double)Rate/(double)Nmeas);
    fflush(stdout);
    parameterfile = fopen(parameterfilename, "a");
    fprintf(parameterfile, "Acceptance Rate was: %e Percent\n", 100.*(double)Rate/(double)Nmeas);
    fclose(parameterfile);
  }

#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_tmp();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  free_monomials();
  if(g_running_phmc) {
    free_bispinor_field();
    free_chi_up_spinor_field();
    free_chi_dn_spinor_field();
  }
  /* End IF PHMC */
/*   if(use_stout_flag == 1) */
/*     free_stout_smear_vars(); */

  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}

static char const rcsid[] = "$Id$";
