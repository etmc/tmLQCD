/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
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
#include "global.h"
#include "start.h"
#include "ranlxs.h"
#include "su3spinor.h"
#include "source_generation.h"
#include "operator.h"
#include "invert_eo.h"
#include "solver/solver.h"
#include "geometry_eo.h"
#include "linalg/convert_eo_to_lexic.h"
#include "measurements.h"
#include "correlators.h"
#include "gettime.h"
#include "DDalphaAMG_interface.h"
#include "read_input.h"
#include "init/init_gauge_tmp.h"

/******************************************************
 *
 * This routine computes the correlators
 * <PP>, <PA> and <PV> (<source sink>)
 * using a stochastic time slice source
 * and only one inversion (actually A_0)
 * 
 * for <AP> we would need another inversion
 *
 *
 *
 ******************************************************/

void correlators_measurement(const int traj, const int id, const int ieo) {
  int i, j, t, tt, t0;
  double *Cpp = NULL, *Cpa = NULL, *Cp4 = NULL;
  double res = 0., respa = 0., resp4 = 0.;
  double atime, etime;
  float tmp;
  operator * optr;
#ifdef TM_USE_MPI
  double mpi_res = 0., mpi_respa = 0., mpi_resp4 = 0.;
  // send buffer for MPI_Gather
  double *sCpp = NULL, *sCpa = NULL, *sCp4 = NULL;
#endif
  FILE *ofs;
  char *filename;
  char *filename_tmp;
  char buf[100], buf2[100];
  spinor phi;
  filename=buf;
  filename_tmp = buf2;

  init_operators();
  if(no_operators < 1) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Warning! no operators defined in input file, cannot perform online correlator mesurements!\n");
    }
    return;
  }
  if(no_operators > 1 && g_proc_id == 0) {
    fprintf(stderr, "Warning! number of operators defined larger than 1, using only the first!\n");
  }
  optr = &operator_list[0];
  // we don't want to do inversion twice for this purpose here
  optr->DownProp = 0;
  if(optr->type != TMWILSON && optr->type != WILSON && optr->type != CLOVER) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Warning! correlator online measurement currently only implemented for TMWILSON, WILSON and CLOVER\n");
      fprintf(stderr, "Cannot perform correlator online measurement!\n");
    }
    return;
  }
  
  if(ranlxs_init == 0) {
    rlxs_init(1, 123456);
  }

  // there are three modes of operation
  // 1) one single time-slice source (default)
  // 2) no_samples time-slice sources on random time-slices
  // 3) one sample on all time-slices
  int max_samples = measurement_list[id].all_time_slices ? 1 : measurement_list[id].no_samples;
  int max_time_slices = measurement_list[id].all_time_slices ? measurement_list[id].max_source_slice : 1;
  for(int sample = 0; sample < max_samples; sample++ ){
    for(int ts = 0; ts < max_time_slices; ts++){

      if( max_samples == 1 && max_time_slices == 1 ){
        sprintf(filename,"%s%06d", "onlinemeas.", traj);
      } else if ( max_samples == 1 && max_time_slices > 1){
        sprintf(filename,"%s%06d.t%03d", "onlinemeas.", traj, ts );
      } else {
        sprintf(filename,"%s%06d.s%03d", "onlinemeas.", traj, sample);
      }
      /* generate random timeslice */
      t0 = ts;
      if( !measurement_list[id].all_time_slices ){
        ranlxs(&tmp, 1);
        t0 = (int)(measurement_list[id].max_source_slice*tmp);
      }
#ifdef TM_USE_MPI
      MPI_Bcast(&t0, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if(g_debug_level > 1 && g_proc_id == 0) {
        printf("# timeslice set to %d (T=%d) for online measurement\n", t0, g_nproc_t*T);
        printf("# online measurements parameters: kappa = %.12f, mu = %.12f\n", optr->kappa, optr->mu/2./optr->kappa);
      }
      atime = gettime();

      int runs = 1;
      if (restoresu3_flag) runs = 2;

#ifdef TM_USE_MPI
      sCpp = (double*) calloc(T, sizeof(double));
      sCpa = (double*) calloc(T, sizeof(double));
      sCp4 = (double*) calloc(T, sizeof(double));
      if(g_mpi_time_rank == 0) {
        Cpp = (double*) calloc(g_nproc_t*T*runs, sizeof(double));
        Cpa = (double*) calloc(g_nproc_t*T*runs, sizeof(double));
        Cp4 = (double*) calloc(g_nproc_t*T*runs, sizeof(double));
      }
#else
      Cpp = (double*) calloc(T*runs, sizeof(double));
      Cpa = (double*) calloc(T*runs, sizeof(double));
      Cp4 = (double*) calloc(T*runs, sizeof(double));
#endif
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
	    		      t0, sample, traj, measurement_list[id].seed);
      optr->sr0 = g_spinor_field[0];
      optr->sr1 = g_spinor_field[1];
      optr->prop0 = g_spinor_field[2];
      optr->prop1 = g_spinor_field[3];

      for( int r = 0; r<runs; r++) {
        
        if (restoresu3_flag) {
          for(int ix=0;ix<VOLUME;ix++) {
            for(int mu=0;mu<4;mu++){
              su3 *v, *w;
              v=&(g_gauge_field[ix][mu]);
              w=&(gauge_tmp[ix][mu]);
              if(r == 0){
                _su3_assign(*v,*w);
              } else {
                restoresu3_in_place(v);
              }
            }
          }
#ifdef TM_USE_MPI
          xchange_gauge(g_gauge_field);
#endif
          mg_update_gauge = 1;
        }

        // op_id = 0, index_start = 0, write_prop = 0
        optr->inverter(0, 0, 0);

        /* now we bring it to normal format */
        /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
        convert_eo_to_lexic(g_spinor_field[DUM_MATRIX], g_spinor_field[2], g_spinor_field[3]);
      
        /* now we sum only over local space for every t */
        for(t = 0; t < T; t++) {
          j = g_ipt[t][0][0][0];
          res = 0.;
          respa = 0.;
          resp4 = 0.;
          for(i = j; i < j+LX*LY*LZ; i++) {
            res += _spinor_prod_re(g_spinor_field[DUM_MATRIX][i], g_spinor_field[DUM_MATRIX][i]);
            _gamma0(phi, g_spinor_field[DUM_MATRIX][i]);
            respa += _spinor_prod_re(g_spinor_field[DUM_MATRIX][i], phi);
            _gamma5(phi, phi);
            resp4 += _spinor_prod_im(g_spinor_field[DUM_MATRIX][i], phi);
          }
          
#if defined TM_USE_MPI
          MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
          res = mpi_res;
          MPI_Reduce(&respa, &mpi_respa, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
          respa = mpi_respa;
          MPI_Reduce(&resp4, &mpi_resp4, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
          resp4 = mpi_resp4;
          sCpp[t] = +res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
          sCpa[t] = -respa/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
          sCp4[t] = +resp4/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
#else
          Cpp[t+g_nproc_t*T*r] = +res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
          Cpa[t+g_nproc_t*T*r] = -respa/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
          Cp4[t+g_nproc_t*T*r] = +resp4/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
#endif
        }
        
#ifdef TM_USE_MPI
        /* some gymnastics needed in case of parallelisation */
        if(g_mpi_time_rank == 0) {
          MPI_Gather(sCpp, T, MPI_DOUBLE, Cpp+g_nproc_t*T*r, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
          MPI_Gather(sCpa, T, MPI_DOUBLE, Cpa+g_nproc_t*T*r, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
          MPI_Gather(sCp4, T, MPI_DOUBLE, Cp4+g_nproc_t*T*r, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        }
#endif
        
        /* and write everything into a file */
        if(g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
          if(runs > 1) {
            sprintf(filename_tmp,"%s.r%02d", filename, r);
            ofs = fopen(filename_tmp, "w");
          } else {
            ofs = fopen(filename, "w");
          }
          fprintf( ofs, "1  1  0  %e  %e\n", Cpp[t0+g_nproc_t*T*r], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "1  1  %d  %e  ", t, Cpp[tt+g_nproc_t*T*r]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cpp[tt+g_nproc_t*T*r]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "1  1  %d  %e  %e\n", t, Cpp[tt+g_nproc_t*T*r], 0.);
          
          fprintf( ofs, "2  1  0  %e  %e\n", Cpa[t0+g_nproc_t*T*r], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "2  1  %d  %e  ", t, Cpa[tt+g_nproc_t*T*r]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cpa[tt+g_nproc_t*T*r]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "2  1  %d  %e  %e\n", t, Cpa[tt+g_nproc_t*T*r], 0.);
          
          fprintf( ofs, "6  1  0  %e  %e\n", Cp4[t0+g_nproc_t*T*r], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "6  1  %d  %e  ", t, Cp4[tt+g_nproc_t*T*r]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cp4[tt+g_nproc_t*T*r]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "6  1  %d  %e  %e\n", t, Cp4[tt+g_nproc_t*T*r], 0.);
          fclose(ofs);
        }
      }
      if(g_mpi_time_rank == 0 && g_proc_coords[0] == 0 && runs == 2) {
        sprintf(filename_tmp,"%s.diff", filename);
        ofs = fopen(filename_tmp, "w");
        fprintf( ofs, "1  1  0  %e  %e\n", Cpp[t0] - Cpp[t0+g_nproc_t*T], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "1  1  %d  %e  ", t, Cpp[tt] - Cpp[tt+g_nproc_t*T]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cpp[tt] - Cpp[tt+g_nproc_t*T]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "1  1  %d  %e  %e\n", t, Cpp[tt] - Cpp[tt+g_nproc_t*T], 0.);
          
          fprintf( ofs, "2  1  0  %e  %e\n", Cpa[t0] - Cpa[t0+g_nproc_t*T], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "2  1  %d  %e  ", t, Cpa[tt] - Cpa[tt+g_nproc_t*T]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cpa[tt] - Cpa[tt+g_nproc_t*T]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "2  1  %d  %e  %e\n", t, Cpa[tt] - Cpa[tt+g_nproc_t*T], 0.);
          
          fprintf( ofs, "6  1  0  %e  %e\n", Cp4[t0] - Cp4[t0+g_nproc_t*T], 0.);
          for(t = 1; t < g_nproc_t*T/2; t++) {
            tt = (t0+t)%(g_nproc_t*T);
            fprintf( ofs, "6  1  %d  %e  ", t, Cp4[tt] - Cp4[tt+g_nproc_t*T]);
            tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
            fprintf( ofs, "%e\n", Cp4[tt] - Cp4[tt+g_nproc_t*T]);
          }
          tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
          fprintf( ofs, "6  1  %d  %e  %e\n", t, Cp4[tt] - Cp4[tt+g_nproc_t*T], 0.);
          fclose(ofs);
        }

#ifdef TM_USE_MPI
      if(g_mpi_time_rank == 0) {
        free(Cpp); free(Cpa); free(Cp4);
      }
      free(sCpp); free(sCpa); free(sCp4);
#else
      free(Cpp); free(Cpa); free(Cp4);
#endif
    }
  } 
  etime = gettime();
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("ONLINE: measurement done int t/s = %1.4e\n", etime - atime);
  }
  return;
}
