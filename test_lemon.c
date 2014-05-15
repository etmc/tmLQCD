/***********************************************************************
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
/*******************************************************************************
*
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
*
*
*******************************************************************************/

#include <lime.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef MPI
# include <mpi.h>
#endif
#include "su3.h"
#include "su3adj.h"
#include <io/params.h>
#include <io/gauge.h>

#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "boundary.h"
#include "global.h"
#include "xchange/xchange.h"
#include "init/init.h"
#include "measure_gauge_action.h"
#include "mpi_init.h"


int main(int argc,char *argv[]) {

  double plaquette_energy;
  paramsXlfInfo *xlfInfo;
  

#ifdef MPI
  
  MPI_Init(&argc, &argv);
#endif
  g_rgi_C1 = 1.; 
  
  /* Read the input file */
  read_input("benchmark.input");
  
  tmlqcd_mpi_init(argc, argv);
  
  
#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

  if(g_proc_id == 0) {
    fprintf(stdout,"The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    printf("# Testing IO routines for gauge-fields\n");
    fflush(stdout);
  }
  
  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

  /* generate a random gauge field */
  start_ranlux(1, 123456);
  random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);

#ifdef MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif

  plaquette_energy = measure_plaquette(g_gauge_field) / (6.*VOLUME*g_nproc);

  if(g_proc_id == 0) {
    printf("# the first plaquette value is %e\n", plaquette_energy);
    printf("# writing with lime first to conf.lime\n");
  }

  /* write with lime first */
  xlfInfo = construct_paramsXlfInfo(plaquette_energy, 0);
  write_lime_gauge_field( "conf.lime", 64, xlfInfo);

#ifdef HAVE_LIBLEMON
  if(g_proc_id == 0) {
    printf("Now we do write with lemon to conf.lemon...\n");
  }
  write_lemon_gauge_field_parallel( "conf.lemon", 64, xlfInfo);


  if(g_proc_id == 0) {
    printf("# now we read with lemon from conf.lime\n");
  }
  read_lemon_gauge_field_parallel("conf.lime", NULL, NULL, NULL);
  plaquette_energy = measure_plaquette(g_gauge_field) / (6.*VOLUME*g_nproc);
  if(g_proc_id == 0) {
    printf("# the plaquette value after lemon read of conf.lime is %e\n", plaquette_energy);
  }

  if(g_proc_id == 0) {
    printf("# now we read with lemon from conf.lemon\n");
  }
  read_lemon_gauge_field_parallel("conf.lemon", NULL, NULL, NULL);
  plaquette_energy = measure_plaquette(g_gauge_field) / (6.*VOLUME*g_nproc);
  if(g_proc_id == 0) {
    printf("# the plaquette value after lemon read of conf.lemon is %e\n", plaquette_energy);
  }

  if(g_proc_id == 0) {
    printf("# now we read with lime from conf.lemon\n");
  }
  read_lime_gauge_field("conf.lemon");
  plaquette_energy = measure_plaquette(g_gauge_field) / (6.*VOLUME*g_nproc);
  if(g_proc_id == 0) {
    printf("# the plaquette value after lime read of conf.lemon is %e\n", plaquette_energy);
  }

  free(xlfInfo);
  if(g_proc_id==0) {
    printf("done ...\n");
  }
#endif

  if(g_proc_id == 0) {
    printf("# now we read with lime from conf.lime\n");
  }
  read_lime_gauge_field("conf.lime", NULL, NULL, NULL);
  plaquette_energy = measure_plaquette(g_gauge_field) / (6.*VOLUME*g_nproc);
  if(g_proc_id == 0) {
    printf("# the plaquette value after lime read of conf.lime is %e\n", plaquette_energy);
  }


#ifdef MPI
  MPI_Finalize();
#endif
  free_gauge_field();
  free_geometry_indices();
  return(0);
}
