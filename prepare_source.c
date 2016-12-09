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
#ifdef TM_USE_MPI
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
#include "ranlxs.h"
#include "su3.h"
#include "operator.h"
#include "linalg_eo.h"
#include "operator/tm_operators_nd.h"
#include "source_generation.h"
#include "prepare_source.h"

void prepare_source(const int nstore, const int isample, const int ix, const int op_id, 
                    const int read_source_flag,
                    const int source_location, const unsigned int seed) {

  FILE * ifs = NULL;
  int is = ix / 3, ic = ix %3, err = 0, rstat=0, t = 0;
  operator * optr = &operator_list[op_id];
  char source_filename[400];
  int source_type = SourceInfo.type;
  float u;
  SourceInfo.nstore = nstore;
  SourceInfo.sample = isample;
  SourceInfo.ix = ix;

  if(optr->type != DBTMWILSON && optr->type != DBCLOVER) {
    SourceInfo.no_flavours = 1;
    /* no volume sources */
    if(source_type == SRC_TYPE_POINT || source_type == SRC_TYPE_TS) {
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
                if(T_global > 99) sprintf(source_filename, "%s.%.4d.%.3d.%.2d", SourceInfo.basename, nstore, t, ix);
                else sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, t, ix);
                if( (ifs = fopen(source_filename, "r")) != NULL) {
                  fclose(ifs);
                  break;
                }
              }
            }
#ifdef TM_USE_MPI
            MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
            SourceInfo.t = t;
          }
          if(T_global > 99) sprintf(source_filename, "%s.%.4d.%.3d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
          else sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
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
        if(T_global > 99) sprintf(source_filename, "%s.%.4d.%.3d.%.2d.inverted", PropInfo.basename, nstore, SourceInfo.t, ix);
        else sprintf(source_filename, "%s.%.4d.%.2d.%.2d.inverted", PropInfo.basename, nstore, SourceInfo.t, ix);
      }
      else {
        if(T_global > 99) sprintf(source_filename, "%s.%.4d.%.3d.inverted", PropInfo.basename, nstore, SourceInfo.t);
        else sprintf(source_filename, "%s.%.4d.%.2d.inverted", PropInfo.basename, nstore, SourceInfo.t);
      }
    }
    else if(source_type == SRC_TYPE_VOL) {
      /* Volume sources */
      if(read_source_flag == 0 || read_source_flag == 2) {
        if(g_proc_id == 0 && g_debug_level > 0) {
          printf("# Preparing 1 flavour volume source\n");
        }
        gaussian_volume_source(g_spinor_field[0], g_spinor_field[1], isample, nstore, 0);
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
    else if(source_type == SRC_TYPE_PION_TS) {
      // If a pion timeslice source has already been inverted for the current sample and gauge configuration,
      // we would like to re-use the same timeslice, which we ensure with the loop below. The reason for doing 
      // this is that we cannot guarantee that the call to ranlxs below is reproducible.
      // Note: source_generation_pion_only reinitialises the RNG with a systematically chosen seed and thus does
      // not suffer from this problem when called below.
      if(SourceInfo.automaticTS) {
        int found = 0;
        if(g_proc_id == 0 && !PropInfo.splitted) {
          for(t = 0; t < g_nproc_t*T; t++) {
            sprintf(source_filename, "%s.%.4d.%.5d.%.2d.inverted", SourceInfo.basename, nstore, isample, t);
            if( (ifs = fopen(source_filename, "r")) != NULL) {
              fclose(ifs);
              found = 1;
              break;
            }
          }
        }
        // chose timeslice randomly
        if(PropInfo.splitted || !found) {
          ranlxs(&u, 1);
          t = (int)(u*g_nproc_t*T);
        }
#ifdef TM_USE_MPI
        MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
        SourceInfo.t = t;
      }
      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("# Preparing 1 flavour Pion TimeSlice at t = %d source\n", SourceInfo.t);
      }
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], SourceInfo.t, isample, nstore, seed);
      sprintf(source_filename, "%s.%.4d.%.5d.%.2d.inverted", PropInfo.basename, nstore, isample, SourceInfo.t);
    }
    else if(source_type == SRC_TYPE_GEN_PION_TS) {
      // Generalised Pion full time slice sources
      if(SourceInfo.automaticTS) {
        // automatic timeslice detection based on an existing forward propagator
        if(g_proc_id == 0) {
          for(t = 0; t < g_nproc_t*T; t++) {
            sprintf(source_filename, "%s.%.4d.%.5d.%.2d.inverted", SourceInfo.basename, nstore, isample, t);
            if( (ifs = fopen(source_filename, "r")) != NULL) {
              fclose(ifs);
              break;
            }
          }
        }
#ifdef TM_USE_MPI
        MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
        SourceInfo.t = t;
      }

      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("# Preparing 1 flavour Generalised Pion TimeSlice at T/2 + t = %d source\n", SourceInfo.t+(g_nproc_t*T)/2);
      }

      sprintf(source_filename, "%s.%.4d.%.5d.%.2d.inverted", SourceInfo.basename, nstore, isample, SourceInfo.t);
      rstat = read_spinor(g_spinor_field[2], g_spinor_field[3], source_filename, 0);
      if(rstat) {
        fprintf(stderr, "Error reading file %s in prepare_source.c.\nUnable to proceed, aborting....\n", source_filename);
        exit(-1);
      }
      extended_pion_source(g_spinor_field[0], g_spinor_field[1], g_spinor_field[2], g_spinor_field[3], 
                           SourceInfo.t, (g_nproc_t*T)/2, 0., 0., 0.);
      sprintf(source_filename, "%s.%.4d.%.5d.%.2d.inverted", PropInfo.basename, nstore, isample, SourceInfo.t);
      // if the generalised pion propagator is to be written to the same file as the source, splitting must be disabled
      if( strcmp(PropInfo.basename,SourceInfo.basename) == 0 ) PropInfo.splitted = 0;
    }
    else { 
      fprintf(stderr, "# source type %d not implemented yet.\nCannot proceed, aborting...\n", source_type);
    }
    optr->sr0 = g_spinor_field[0];
    optr->sr1 = g_spinor_field[1];
    optr->prop0 = g_spinor_field[2];
    optr->prop1 = g_spinor_field[3];


    /* If the solver is _not_ CG we might read in */
    /* here some better guess                     */
    /* This also works for re-iteration           */
    if (optr->solver != CG && optr->solver != PCG && optr->solver != MIXEDCG && optr->solver != RGMIXEDCG) {
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
          mul_r(optr->prop1, 1. / (2*optr->kappa), optr->prop1, VOLUME / 2);
          mul_r(optr->prop0, 1. / (2*optr->kappa), optr->prop0, VOLUME / 2);
        }

        if (err != 0) {
          zero_spinor_field(optr->prop0, VOLUME / 2);
          zero_spinor_field(optr->prop1, VOLUME / 2);
        }
      }
      else {
        zero_spinor_field(optr->prop0, VOLUME / 2);
        zero_spinor_field(optr->prop1, VOLUME / 2);
      }
    }
    else {
      zero_spinor_field(optr->prop0, VOLUME / 2);
      zero_spinor_field(optr->prop1, VOLUME / 2);
    }
    /*     if(optr->even_odd_flag) { */
    /*       assign(optr->sr0, g_spinor_field[0], VOLUME/2); */
    /*       assign(optr->sr1, g_spinor_field[1], VOLUME/2); */
    /*     } */
    /*     else { */
    /*       convert_eo_to_lexic(optr->sr0, g_spinor_field[0], g_spinor_field[1]); */
    /*     } */
  }
  else { /* for the ND 2 flavour twisted operator */
    SourceInfo.no_flavours = 2;
    zero_spinor_field(g_spinor_field[0], VOLUME/2);
    zero_spinor_field(g_spinor_field[1], VOLUME/2);
    if(source_type == SRC_TYPE_POINT || source_type == SRC_TYPE_TS) {
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
          if(T_global > 99) sprintf(source_filename, "%s.%.4d.%.3d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
          else sprintf(source_filename, "%s.%.4d.%.2d.%.2d", SourceInfo.basename, nstore, SourceInfo.t, ix);
        }
        else {
          sprintf(source_filename,"%s", SourceInfo.basename);
        }
        if(g_proc_id == 0) {
          printf("# Trying to read source from %s\n", source_filename);
        }
        if(read_spinor(g_spinor_field[2], g_spinor_field[3], source_filename, 0) != 0) {
          fprintf(stderr, "Error reading source! Aborting...\n");
#ifdef TM_USE_MPI
          MPI_Abort(MPI_COMM_WORLD, 1);
          MPI_Finalize();
#endif
          exit(-1);
        }
      }
    }
    else if(source_type == SRC_TYPE_VOL) {
      /* Volume sources */
      if(g_proc_id == 0 && g_debug_level > 0) {
        printf("# Preparing 2 flavour volume source\n");
      }
      gaussian_volume_source(g_spinor_field[0], g_spinor_field[1],
                             isample, nstore, 1);
      gaussian_volume_source(g_spinor_field[2], g_spinor_field[3],
                             isample, nstore, 2);
    }
    mul_one_pm_itau2(g_spinor_field[4], g_spinor_field[6], g_spinor_field[0], g_spinor_field[2], +1., VOLUME/2);
    mul_one_pm_itau2(g_spinor_field[5], g_spinor_field[7], g_spinor_field[1], g_spinor_field[3], +1., VOLUME/2);
    assign(g_spinor_field[0], g_spinor_field[4], VOLUME/2);
    assign(g_spinor_field[1], g_spinor_field[5], VOLUME/2);
    assign(g_spinor_field[2], g_spinor_field[6], VOLUME/2);
    assign(g_spinor_field[3], g_spinor_field[7], VOLUME/2);

    optr->sr0 = g_spinor_field[0];
    optr->sr1 = g_spinor_field[1];
    optr->sr2 = g_spinor_field[2];
    optr->sr3 = g_spinor_field[3];
    optr->prop0 = g_spinor_field[4];
    optr->prop1 = g_spinor_field[5];
    optr->prop2 = g_spinor_field[6];
    optr->prop3 = g_spinor_field[7];
  }
  return;
}
