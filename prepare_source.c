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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver/solver.h"
#include "start.h"
#include "ranlxd.h"
#include "su3.h"
#include "operator.h"
#include "linalg_eo.h"
#include "operator/tm_operators_nd.h"
#include "source_generation.h"
#include "prepare_source.h"

void prepare_source(const int nstore, const int isample, const int ix, const int op_id, 
                    const int read_source_flag,
                    const int source_location) {

  FILE * ifs = NULL;
  int is = ix / 3, ic = ix %3, err = 0, rstat=0, t = 0;
  operator * optr = &operator_list[op_id];
  char source_filename[100];
  int source_type = SourceInfo.type;
  static int nstore_ = -1;
  static int isample_ = -1;
  static int ix_ = -1;
  static int op_id_ = -1;
  int Vol;
  
  if(optr->generic_noeo){
    Vol = VOLUME;    
  }
  else{
    Vol = VOLUME / 2;
  }
  
  SourceInfo.nstore = nstore;
  SourceInfo.sample = isample;
  SourceInfo.ix = ix;

  if(optr->type != DBTMWILSON && optr->type != DBCLOVER) {
    SourceInfo.no_flavours = 1;
    /* no volume sources */
    if(source_type != 1) {
      /* either "Don't read inversion source from file" or                    */
      /* "Don't read inversion source from file, but save the one generated" */
      if (read_source_flag == 0 || read_source_flag == 2) {
        if (source_location == 0) {
          source_spinor_field(g_spinor_field[0], g_spinor_field[1], is, ic);
        }
        else {
          source_spinor_field_point_from_file(g_spinor_field[0], g_spinor_field[1], is, ic, source_location);
        }
      }
      /* "Read inversion source from file" */
      else {
        if (SourceInfo.splitted) {
	  /* timeslice needs to be put into filename */
	  if(SourceInfo.automaticTS) {
	    /* automatic timeslice detection */
	    if(g_proc_id == 0) {
	      for(t = 0; t < g_nproc_t*T; t++) {
		sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, t, ix);
		if( (ifs = fopen(source_filename, "r")) != NULL) {
		  fclose(ifs);
		  break;
		}
	      }
	    }
#ifdef MPI
	    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
	    SourceInfo.t = t;
	  }
          sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
          if (g_cart_id == 0) {
            printf("# Trying to read source from %s\n", source_filename);
          }
          rstat = read_spinor(g_spinor_field[0], g_spinor_field[1], source_filename, 0);
        }
        else {
          sprintf(source_filename, "%s", SourceInfo.basename);
          if (g_cart_id == 0) {
            printf("# Trying to read source no %d from %s\n", ix, source_filename);
          }
          rstat = read_spinor(g_spinor_field[0], g_spinor_field[1], source_filename, ix);
        }
        if(rstat) {
          fprintf(stderr, "Error reading file %s in prepare_source.c\nUnable to proceed, aborting....\n", source_filename);
          exit(-1);
        }
      }
      if (PropInfo.splitted) {
        sprintf(source_filename, "%s.%.4d.%.2d.%.2d.inverted", PropInfo.basename, nstore, SourceInfo.t, ix);
      }
      else {
        sprintf(source_filename, "%s.%.4d.%.2d.inverted", PropInfo.basename, nstore, SourceInfo.t);
      }
    }
    else if(source_type == 1) {
      /* Volume sources */
      if(read_source_flag == 0 || read_source_flag == 2) {

        if(optr->generic_noeo){
	  if(g_proc_id == 0 && g_debug_level > 0) {
	    printf("# Preparing 1 flavour no-EO volume source\n");
	  }	  
          gaussian_volume_source_noeo(g_spinor_field[0], isample, nstore, 0);
	}
	else{
	  if(g_proc_id == 0 && g_debug_level > 0) {
	    printf("# Preparing 1 flavour volume source\n");
	  }	  
          gaussian_volume_source(g_spinor_field[0], g_spinor_field[1], isample, nstore, 0);
	}
      }
      else {
        sprintf(source_filename, "%s.%.4d.%.5d", SourceInfo.basename, nstore, isample);
        if (g_cart_id == 0) {
          printf("# Trying to read source from %s\n", source_filename);
        }
        rstat = read_spinor(g_spinor_field[0], g_spinor_field[1], source_filename, 0);
        if(rstat) {
          fprintf(stderr, "Error reading file %s in prepare_source.c.\nUnable to proceed, aborting....\n", source_filename);
          exit(-1);
        }
      }
      sprintf(source_filename, "%s.%.4d.%.5d.inverted", PropInfo.basename, nstore, isample);
    }
    
    optr->sr0 = g_spinor_field[0];
    optr->sr1 = g_spinor_field[1];
    optr->prop0 = g_spinor_field[2];
    optr->prop1 = g_spinor_field[3];  

    /* If the solver is _not_ CG we might read in */
    /* here some better guess                     */
    /* This also works for re-iteration           */
    if (optr->solver != CG && optr->solver != PCG) {
      ifs = fopen(source_filename, "r");
      if (ifs != NULL) {
        if (g_cart_id == 0) {
          printf("# Trying to read guess from file %s\n", source_filename);
          fflush(stdout);
        }
        fclose(ifs);
        err = 0;
        /* iter = get_propagator_type(source_filename); */
        rstat = read_spinor(optr->prop0, optr->prop1, source_filename, (PropInfo.splitted ? 0 : ix));
        if(rstat) {
          fprintf(stderr, "Error reading file %s in prepare_source.c, rstat = %d\n", source_filename, rstat);
          exit(-1);
        }
        if (g_kappa != 0.) {
          mul_r(optr->prop1, 1. / (2*optr->kappa), optr->prop1, Vol);
          mul_r(optr->prop0, 1. / (2*optr->kappa), optr->prop0, Vol);
        }

        if (err != 0) {
          zero_spinor_field(optr->prop0, Vol);
          zero_spinor_field(optr->prop1, Vol);
        }
      }
      else {
        zero_spinor_field(optr->prop0, Vol);
        zero_spinor_field(optr->prop1, Vol);
      }
    }
    else {
        zero_spinor_field(optr->prop0, Vol);
        zero_spinor_field(optr->prop1, Vol);
    }
    /*     if(optr->even_odd_flag) { */
    /*       assign(optr->sr0, g_spinor_field[0], VOLUME/2); */
    /*       assign(optr->sr1, g_spinor_field[1], VOLUME/2); */
    /*     } */
    /*     else { */
    /*       convert_eo_to_lexic(optr->sr0, g_spinor_field[0], g_spinor_field[1]); */
    /*     } */
    
    if(optr->generic_noeo){
      /* We need to put the second source/propagator field
	to NULL in order to use the correct (i.e. the lexicographic)
	binary write functions in op_invert FIXME
      */
      optr->sr1 = NULL;
      optr->prop1 = NULL;
    }
  }
  else { /* for the ND 2 flavour twisted operator */
    SourceInfo.no_flavours = 2;
    zero_spinor_field(g_spinor_field[0], Vol);
    zero_spinor_field(g_spinor_field[1], Vol);
    if(source_type != 1) {
      if(read_source_flag == 0 || read_source_flag == 2) {
        if(source_location == 0) {
          source_spinor_field(g_spinor_field[2], g_spinor_field[3], is, ic);
        }
        else {
          source_spinor_field_point_from_file(g_spinor_field[2], g_spinor_field[3], 
                              is, ic, source_location);
        }
          }
          else {
        if(SourceInfo.splitted) {
          sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
        }
        else {
          sprintf(source_filename,"%s", SourceInfo.basename);
        }
        if(g_proc_id == 0) {
          printf("# Trying to read source from %s\n", source_filename);
        }
        if(read_spinor(g_spinor_field[2], g_spinor_field[3], source_filename, 0) != 0) {
          fprintf(stderr, "Error reading source! Aborting...\n");
#ifdef MPI
          MPI_Abort(MPI_COMM_WORLD, 1);
          MPI_Finalize();
#endif
          exit(-1);
        }
      }
    }
    else if(source_type == 1) {
      /* Volume sources */
      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("# Preparing 2 flavour volume source\n");
      }
      if(optr->generic_noeo){
	gaussian_volume_source_noeo(g_spinor_field[0], 
			      isample, nstore, 1);
	gaussian_volume_source_noeo(g_spinor_field[2], 
			      isample, nstore, 2);	
      }
      else{
	gaussian_volume_source(g_spinor_field[0], g_spinor_field[1],
			      isample, nstore, 1);
	gaussian_volume_source(g_spinor_field[2], g_spinor_field[3],
			      isample, nstore, 2);
      }
    }
    mul_one_pm_itau2(g_spinor_field[4], g_spinor_field[6], g_spinor_field[0], g_spinor_field[2], +1., Vol);
    mul_one_pm_itau2(g_spinor_field[5], g_spinor_field[7], g_spinor_field[1], g_spinor_field[3], +1., Vol);
    assign(g_spinor_field[0], g_spinor_field[4], Vol);
    assign(g_spinor_field[1], g_spinor_field[5], Vol);
    assign(g_spinor_field[2], g_spinor_field[6], Vol);
    assign(g_spinor_field[3], g_spinor_field[7], Vol);

    optr->sr0 = g_spinor_field[0];
    optr->sr1 = g_spinor_field[1];
    optr->sr2 = g_spinor_field[2];
    optr->sr3 = g_spinor_field[3];
    optr->prop0 = g_spinor_field[4];
    optr->prop1 = g_spinor_field[5];
    optr->prop2 = g_spinor_field[6];
    optr->prop3 = g_spinor_field[7];
    if(optr->generic_noeo){
      /* We need to put the second source/propagator field
	to NULL in order to use the correct (i.e. the lexicographic)
	binary write functions in op_invert FIXME
      */
      optr->sr1 = NULL;
      optr->sr3 = NULL;      
      optr->prop1 = NULL;
      optr->prop3 = NULL;      
    }   
  }
  nstore_ = nstore;
  isample_ = isample;
  ix_ = ix;
  op_id_ = op_id;
  return;
}
