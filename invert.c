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
 * invert for twisted mass QCD
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include"lime.h"
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "getopt.h"
#include "linalg_eo.h"
#include "geometry_eo.h"
#include "start.h"
/*#include "eigenvalues.h"*/
#include "observables.h"
#ifdef MPI
#include "xchange.h"
#endif
#include <io/utils.h>
#include "read_input.h"
#include "mpi_init.h"
#include "sighandler.h"
#include "boundary.h"
#include "solver/solver.h"
#include "init_gauge_field.h"
#include "init_geometry_indices.h"
#include "init_spinor_field.h"
#include "init_moment_field.h"
#include "init_dirac_halfspinor.h"
#include "init_bispinor_field.h"
#include "init_chi_spinor_field.h"
#include "xchange_halffield.h"
#include "stout_smear.h"
#include "invert_eo.h"
#include "monomial.h"
#include "ranlxd.h"
#include "phmc.h"
#include "D_psi.h"
#include "little_D.h"
#include "reweighting_factor.h"
#include "linalg/convert_eo_to_lexic.h"
#include "block.h"
#include "operator.h"
#include "sighandler.h"
#include "solver/dfl_projector.h"
#include "solver/generate_dfl_subspace.h"
#include "prepare_source.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>

void usage()
{
  fprintf(stdout, "Inversion for EO preconditioned Wilson twisted mass QCD\n");
  fprintf(stdout, "Version %s \n\n", PACKAGE_VERSION);
  fprintf(stdout, "Please send bug reports to %s\n", PACKAGE_BUGREPORT);
  fprintf(stdout, "Usage:   invert [options]\n");
  fprintf(stdout, "Options: [-f input-filename]\n");
  fprintf(stdout, "         [-o output-filename]\n");
  fprintf(stdout, "         [-h|-? this help]\n");
  exit(0);
}

extern int nstore;
int check_geometry();

int main(int argc, char *argv[])
{
  FILE *parameterfile = NULL;
  int c, j, ix = 0, isample = 0, op_id = 0;
  char * filename = NULL;
  char datafilename[50];
  char parameterfilename[50];
  char conf_filename[50];
  char * input_filename = NULL;
  double plaquette_energy;

#ifdef _KOJAK_INST
#pragma pomp inst init
#pragma pomp inst begin(main)
#endif
  

#if (defined SSE || defined SSE2 || SSE3)
  signal(SIGILL, &catch_ill_inst);
#endif

  DUM_DERI = 8;
  /* DUM_DERI + 2 is enough (not 7) */
  DUM_SOLVER = DUM_DERI + 5;
  DUM_MATRIX = DUM_SOLVER + 8;
  /* DUM_MATRIX + 2 is enough (not 6) */
  NO_OF_SPINORFIELDS = DUM_MATRIX + 2;

  verbose = 1;
  g_use_clover_flag = 0;

#ifdef MPI
  MPI_Init(&argc, &argv);
#endif

  while ((c = getopt(argc, argv, "h?f:o:")) != -1) {
    switch (c) {
      case 'f':
        input_filename = calloc(200, sizeof(char));
        strcpy(input_filename, optarg);
        break;
      case 'o':
        filename = calloc(200, sizeof(char));
        strcpy(filename, optarg);
        break;
      case 'h':
      case '?':
      default:
        usage();
        break;
    }
  }
  if (input_filename == NULL) {
    input_filename = "hmc.input";
  }
  if (filename == NULL) {
    filename = "output";
  }

  /* Read the input file */
  if( (j = read_input(input_filename)) != 0) {
    fprintf(stderr, "Could not find input file: %s\nAborting...\n", input_filename);
    exit(-1);
  }

  /* this DBW2 stuff is not needed for the inversion ! */
  if (g_dflgcr_flag == 1) {
    even_odd_flag = 0;
  }
  g_rgi_C1 = 0;
  if (Nsave == 0) {
    Nsave = 1;
  }

  if (g_running_phmc) {
    NO_OF_SPINORFIELDS = DUM_MATRIX + 8;
  }

  tmlqcd_mpi_init(argc, argv);

  g_dbw2rand = 0;

  /* starts the single and double precision random number */
  /* generator                                            */
  start_ranlux(rlxd_level, random_seed);

  /* we need to make sure that we don't have even_odd_flag = 1 */
  /* if any of the operators doesn't use it                    */
  /* in this way even/odd can still be used by other operators */
  for(j = 0; j < no_operators; j++) if(!operator_list[j].even_odd_flag) even_odd_flag = 0;

#ifndef MPI
  g_dbw2rand = 0;
#endif

#ifdef _GAUGE_COPY
  j = init_gauge_field(VOLUMEPLUSRAND, 1);
#else
  j = init_gauge_field(VOLUMEPLUSRAND, 0);
#endif
  if (j != 0) {
    fprintf(stderr, "Not enough memory for gauge_fields! Aborting...\n");
    exit(-1);
  }
  j = init_geometry_indices(VOLUMEPLUSRAND);
  if (j != 0) {
    fprintf(stderr, "Not enough memory for geometry indices! Aborting...\n");
    exit(-1);
  }
  if (no_monomials > 0) {
    if (even_odd_flag) {
      j = init_monomials(VOLUMEPLUSRAND / 2, even_odd_flag);
    }
    else {
      j = init_monomials(VOLUMEPLUSRAND, even_odd_flag);
    }
    if (j != 0) {
      fprintf(stderr, "Not enough memory for monomial pseudo fermion  fields! Aborting...\n");
      exit(0);
    }
  }
  if (even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND / 2, NO_OF_SPINORFIELDS);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, NO_OF_SPINORFIELDS);
  }
  if (j != 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(-1);
  }

  if (g_running_phmc) {
    j = init_chi_spinor_field(VOLUMEPLUSRAND / 2, 20);
    if (j != 0) {
      fprintf(stderr, "Not enough memory for PHMC Chi fields! Aborting...\n");
      exit(0);
    }
  }

  g_mu = g_mu1;
  if (g_cart_id == 0) {
    /*construct the filenames for the observables and the parameters*/
    strcpy(datafilename, filename);
    strcat(datafilename, ".data");
    strcpy(parameterfilename, filename);
    strcat(parameterfilename, ".para");

    parameterfile = fopen(parameterfilename, "w");
    write_first_messages(parameterfile, 1);
    fclose(parameterfile);
  }

  /* define the geometry */
  geometry();

  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  phmc_invmaxev = 1.;

  init_operators();

  /* this could be maybe moved to init_operators */
#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if (j != 0) {
    fprintf(stderr, "Not enough memory for halffield! Aborting...\n");
    exit(-1);
  }
  if (g_sloppy_precision_flag == 1) {
    j = init_dirac_halfspinor32();
    if (j != 0)
    {
      fprintf(stderr, "Not enough memory for 32-Bit halffield! Aborting...\n");
      exit(-1);
    }
  }
#  if (defined _PERSISTENT)
  if (even_odd_flag)
    init_xchange_halffield();
#  endif
#endif

  for (j = 0; j < Nmeas; j++) {
    sprintf(conf_filename, "%s.%.4d", gauge_input_filename, nstore);
    if (g_cart_id == 0) {
      printf("Reading gauge field from file %s\n", conf_filename);
      fflush(stdout);
    }
    if( (j = read_gauge_field(conf_filename)) !=0) {
      fprintf(stderr, "error %d while reading gauge field from %s\n Aborting...\n", j, conf_filename);
      exit(-2);
    }


    if (g_cart_id == 0) {
      printf("done!\n");
      fflush(stdout);
    }
#ifdef MPI
    xchange_gauge();
#endif

    /*compute the energy of the gauge field*/
    plaquette_energy = measure_gauge_action();

    if (g_cart_id == 0) {
      printf("The plaquette value is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
      fflush(stdout);
    }

    if (use_stout_flag == 1) {
      if (stout_smear_gauge_field(stout_rho , stout_no_iter) != 0)
        exit(1) ;

      plaquette_energy = measure_gauge_action();

      if (g_cart_id == 0) {
        printf("The plaquette value after stouting is %e\n", plaquette_energy / (6.*VOLUME*g_nproc));
        fflush(stdout);
      }
    }

    if (reweighting_flag == 1) {
      reweighting_factor(reweighting_samples, nstore);
    }

    /* Compute minimal eigenvalues, if wanted */
    if (compute_evs != 0) {
      eigenvalues(&no_eigenvalues, 5000, eigenvalue_precision,
                  0, compute_evs, nstore, even_odd_flag);
    }
    if (phmc_compute_evs != 0) {
#ifdef MPI
      MPI_Finalize();
#endif
      return(0);
    }

    /* move to operators as well */
    if (g_dflgcr_flag == 1) {
      /* set up deflation blocks */
      init_blocks(1, 1, 2, 1, g_N_s);

      /* the can stay here for now, but later we probably need */
      /* something like init_dfl_solver called somewhere else  */
      /* create set of approximate lowest eigenvectors ("global deflation subspace") */

      /*       g_mu = 0.; */
      /*       boundary(0.125); */
      generate_dfl_subspace(g_N_s, VOLUME);
      /*       boundary(g_kappa); */
      /*       g_mu = g_mu1; */

      /* Compute little Dirac operators */
      /*       alt_block_compute_little_D(); */
      if (g_debug_level > -1)
      {
        check_projectors();
      }
      if (g_debug_level > -1)
      {
        check_little_D_inversion();
      }

    }
    if(SourceInfo.type == 1) {
      index_start = 0;
      index_end = 1;
    }

    for(isample = 0; isample < no_samples; isample++) {
      for (ix = index_start; ix < index_end; ix++) {
	for(op_id = 0; op_id < no_operators; op_id++) {
	  boundary(operator_list[op_id].kappa);
          g_kappa = operator_list[op_id].kappa; 
	  g_mu = 0.;
	  /* we use g_spinor_field[0-7] for sources and props for the moment */
	  /* 0-3 in case of 1 flavour  */
	  /* 0-7 in case of 2 flavours */
	  prepare_source(nstore, isample, ix, op_id, 
			 read_source_flag,
			 source_location);
	  operator_list[op_id].inverter(op_id, index_start);
	}
      }
    }
    nstore += Nsave;
  }
  
#ifdef MPI
  MPI_Finalize();
#endif

  free_blocks();
  free_dfl_subspace();
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
  free_chi_spinor_field();
  return(0);
#ifdef _KOJAK_INST
#pragma pomp inst end(main)
#endif
}
