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
# include<tmlqcd_config.h>
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



#define TM_OMEAS_FILENAME_LENGTH 100

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
void light_correlators_measurement(const int traj, const int id, const int ieo) {
  tm_stopwatch_push(&g_timers, __func__, "");
  int i, j, t, tt, t0;
  double *Cpp = NULL, *Cpa = NULL, *Cp4 = NULL;
  double res = 0., respa = 0., resp4 = 0.;
  // double atime, etime;
  float tmp;
  operator * optr;
#ifdef TM_USE_MPI
  double mpi_res = 0., mpi_respa = 0., mpi_resp4 = 0.;
  // send buffer for MPI_Gather
  double *sCpp = NULL, *sCpa = NULL, *sCp4 = NULL;
#endif
  FILE *ofs;
  char filename[TM_OMEAS_FILENAME_LENGTH];
  spinor phi;

  init_operators();
  if(no_operators < 1) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Warning! no operators defined in input file, cannot perform online correlator mesurements!\n");
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");
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
    tm_stopwatch_pop(&g_timers, 0, 0, "");
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
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, 
                 "%s%06d", "onlinemeas." ,traj);
      } else if ( max_samples == 1 && max_time_slices > 1){
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, 
                 "%s.t%03d.%06d", "onlinemeas", ts, traj );
      } else {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH,
                 "%s.s%03d.%06d", "onlinemeas", sample, traj);
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
      //atime = gettime();

#ifdef TM_USE_MPI
      sCpp = (double*) calloc(T, sizeof(double));
      sCpa = (double*) calloc(T, sizeof(double));
      sCp4 = (double*) calloc(T, sizeof(double));
      if(g_mpi_time_rank == 0) {
        Cpp = (double*) calloc(g_nproc_t*T, sizeof(double));
        Cpa = (double*) calloc(g_nproc_t*T, sizeof(double));
        Cp4 = (double*) calloc(g_nproc_t*T, sizeof(double));
      }
#else
      Cpp = (double*) calloc(T, sizeof(double));
      Cpa = (double*) calloc(T, sizeof(double));
      Cp4 = (double*) calloc(T, sizeof(double));
#endif
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
	    		      t0, sample, traj, measurement_list[id].seed);
      optr->sr0 = g_spinor_field[0];
      optr->sr1 = g_spinor_field[1];
      optr->prop0 = g_spinor_field[2];
      optr->prop1 = g_spinor_field[3];

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
        Cpp[t] = +res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
        Cpa[t] = -respa/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
        Cp4[t] = +resp4/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/2./optr->kappa/optr->kappa;
#endif
      }

#ifdef TM_USE_MPI
      /* some gymnastics needed in case of parallelisation */
      if(g_mpi_time_rank == 0) {
        MPI_Gather(sCpp, T, MPI_DOUBLE, Cpp, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCpa, T, MPI_DOUBLE, Cpa, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCp4, T, MPI_DOUBLE, Cp4, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
      }
#endif

      /* and write everything into a file */
      if(g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
        ofs = fopen(filename, "w");
        fprintf( ofs, "1  1  0  %e  %e\n", Cpp[t0], 0.);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "1  1  %d  %e  ", t, Cpp[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cpp[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "1  1  %d  %e  %e\n", t, Cpp[tt], 0.);

        fprintf( ofs, "2  1  0  %e  %e\n", Cpa[t0], 0.);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "2  1  %d  %e  ", t, Cpa[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cpa[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "2  1  %d  %e  %e\n", t, Cpa[tt], 0.);

        fprintf( ofs, "6  1  0  %e  %e\n", Cp4[t0], 0.);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "6  1  %d  %e  ", t, Cp4[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cp4[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "6  1  %d  %e  %e\n", t, Cp4[tt], 0.);
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
    } // for(max_time_slices)
  } // for(max_samples)
  //etime = gettime();
  //if(g_proc_id == 0 && g_debug_level > 0) {
  //  printf("ONLINE: measurement done int t/s = %1.4e\n", etime - atime);
  //}
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}


// Function to allocate memory for a pointer of pointers for an array of d dimensions
void* callocMultiDimensional(int* sizes, int numDimensions, size_t elementSize) {
    if (numDimensions == 1) {
        void* array = calloc(sizes[0], sizeof(elementSize));
        
        if (array == NULL) {
            printf("Memory allocation failed for dimension 0!");
            return NULL;
        }
        
        return array;
    }
    
    void** pointerArray = (void**)calloc(sizes[0], sizeof(void*));
    
    if (pointerArray == NULL) {
        printf("Memory allocation failed for dimension %d!", numDimensions - 1);
        return NULL;
    }
    
    for (int i = 0; i < sizes[0]; i++) {
        pointerArray[i] = callocMultiDimensional(sizes + 1, numDimensions - 1);
        
        if (pointerArray[i] == NULL) {
            printf("Memory allocation failed for subarray %d in dimension %d!", i, numDimensions - 1);
            // Free memory allocated so far
            for (int j = 0; j < i; j++) {
                free(pointerArray[j]);
            }
            free(pointerArray);
            return NULL;
        }
    }
    
    return pointerArray;
}


/******************************************************
 *
 * This routine computes the correlator matrix <G1,G2> (<source sink>),
 * where G1 and G2 are 2 gamma matrices from the set {1, gamma5} --> {S,P} (scalar and pseudoscalar
 *current), <S,S>, <S,P>, <P,P> for all the choices of h1 and h2 (see eq. 18 and 20 of
 *https://arxiv.org/pdf/1005.2042.pdf):
 *
 * C = Tr[S_l * Gamma_1 * S_h * Gamma_2]
 * where S_l, S_h are the propagators of the light and heavy quark of the meson (e.g. K --> u, )
 *
 * i1, i2 = operator indices from the list of the operators in the input file.
 *
 *
 ******************************************************/
void heavy_correlators_measurement(const int traj, const int id, const int ieo, const int i1,
                                          const int i2) {
  tm_stopwatch_push(&g_timers, __func__, "");  // timer for profiling

  spinor phi;                                  // dummy spinor
  int i, j, t, tt, t0;                         // dummy indices
  float tmp;                                   // dummy variable

  
   /* light-light correlators: dummy variables */

  // array of the values of the correlator C(t)
  double *Cpp = NULL, *Cpa = NULL, *Cp4 = NULL;
  // result of the accumulation of MPI partial sums
  double res = 0., respa = 0., resp4 = 0.;
#ifdef TM_USE_MPI
  double mpi_res = 0., mpi_respa = 0., mpi_resp4 = 0.;
  // send buffer for MPI_Gather
  double *sCpp = NULL, *sCpa = NULL, *sCp4 = NULL;
#endif

#ifdef TM_USE_MPI
  sCpp = (double *)calloc(T, sizeof(double));
  sCpa = (double *)calloc(T, sizeof(double));
  sCp4 = (double *)calloc(T, sizeof(double));
  if (g_mpi_time_rank == 0) {
    Cpp = (double *)calloc(g_nproc_t * T, sizeof(double));
    Cpa = (double *)calloc(g_nproc_t * T, sizeof(double));
    Cp4 = (double *)calloc(g_nproc_t * T, sizeof(double));
  }
#else
  Cpp = (double *)calloc(T, sizeof(double));
  Cpa = (double *)calloc(T, sizeof(double));
  Cp4 = (double *)calloc(T, sizeof(double));
#endif

  /* heavy-light correlators variables */
  double eta_Gamma[2] = {1.0, -1.0}; // sign change compensation from daggering bilinear with Gamma = 1,gamma_5

  // the number of independent correlators is 16 = (2*2)_{h_1 h_2} * (2*2)_{Gamma_1 Gamma_2}
  // hi: c,s 
  // Gamma = 1, gamma_5
  double *****C_hihj_g1g2 = (double *****) callocMultiDimensional({2,2,2,2,T}, 5, sizeof(double));
  double ****res_hihj_g1g2 = (double ****) callocMultiDimensional({2,2,2,2}, 4, sizeof(double));
#ifdef TM_USE_MPI
  double ****mpi_res_hihj_g1g2 = (double ****) callocMultiDimensional({2,2,2,2}, 4, sizeof(double));
  // send buffer for MPI_Gather
  double *****sC_hihj_g1g2 = (double *****) callocMultiDimensional({2,2,2,2,T}, 5, sizeof(double));
#endif

  FILE *ofs;                                // output file stream
  char filename[TM_OMEAS_FILENAME_LENGTH];  // file path

  /* building the stochastic propagator spinors needed for the correlators */

  init_operators();  // initialize the operators in the action (if not already done)
  if (no_operators < 1) {
    if (g_proc_id == 0) {
      // we don't want the simulation to stop, we can do the measurements offline afterwards
      fprintf(stderr,
              "Warning! no operators defined in input file, cannot perform online correlator "
              "mesurements!\n");
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");  // timer for profiling
    return;
  }

  // selecting the operators i1 and i2 from the list of operators initialized before
  operator* optr1 = & operator_list[i1]; // light doublet
  operator* optr2 = & operator_list[i2]; // heavy c,s doublet

  bool b1 = (optr->type != TMWILSON && optr->type != WILSON && optr->type != CLOVER);
  bool b2 = (optr2->type != DBTMWILSON && optr2->type != DBCLOVER);
  if(b1 || b2) {
    if (g_proc_id == 0) {
      if(b1){
        fprintf(stderr, "Warning! optr1 should be one of the following: TMWILSON, WILSON and CLOVER\n");
        fprintf(stderr, "Cannot perform correlator online measurement!\n");
      }
      if (b2) {
        fprintf(stderr, "Warning! optr2 should be one of the following: DBTMWILSON, DBCLOVER\n");
      }
      fprintf(stderr, "Cannot perform correlator online measurement!\n");
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");  // timer for profiling
    return;
  }

  if (ranlxs_init == 0) {
    rlxs_init(1, 123456);  // initializing random number generator RANLUX
  }

  // there are three modes of operation
  // 1) one single time-slice source (default)
  // 2) no_samples time-slice sources on random time-slices
  // 3) one sample on all time-slices
  int max_samples = measurement_list[id].all_time_slices ? 1 : measurement_list[id].no_samples;
  int max_time_slices =
      measurement_list[id].all_time_slices ? measurement_list[id].max_source_slice : 1;
  for (int sample = 0; sample < max_samples; sample++) {
    for (int ts = 0; ts < max_time_slices; ts++) {
      // setting output filename
      if (max_samples == 1 && max_time_slices == 1) {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s%06d", "heavy_mesons.", traj);
      } else if (max_samples == 1 && max_time_slices > 1) {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s.t%03d.%06d", "heavy_mesons", ts, traj);
      } else {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s.s%03d.%06d", "heavy_mesons", sample, traj);
      }

      /* generate random timeslice */
      t0 = ts;
      if (!measurement_list[id].all_time_slices) {
        ranlxs(&tmp, 1);
        t0 = (int)(measurement_list[id].max_source_slice * tmp);
      }
#ifdef TM_USE_MPI
      MPI_Bcast(&t0, 1, MPI_INT, 0, MPI_COMM_WORLD);  // broadcast t0 to all the ranks
#endif
      if (g_debug_level > 1 && g_proc_id == 0) {
        printf("# timeslice set to %d (T=%d) for online measurement\n", t0, g_nproc_t * T);
        printf("# online measurements parameters: kappa = %.12f, mubar = %.12f, epsbar = %.12f\n",
               optr1->kappa, optr1->mubar, optr1->epsbar);
      }

      ///////////// OLD COMMENT
      ///////////// the gamma matrix in front specifies if we are working with 
      ///////////// psi_1 or psi_2 of eqs. (2.7, 2.8) https://arxiv.org/pdf/0804.1501.pdf
      //////////// END OLD COMMENT
      // even-odd spinor fields for the light and heavy doublet correlators
      // INDICES: source+propagator, doublet, spin dilution index, even-odd, flavor, position (+ Dirac, color)
      // psi = (psi[(s,p)][phi][beta][eo][f][x])[alpha][c] // last 3 indices are proper of the spinor struct
      // Note: propagator in the sense that it is D^{-1}*source after the inversion

      spinor ******arr_eo_spinor = (spinor ******)callocMultiDimensional({2,2,4,2,2,VOLUME/2}, 6, sizeof(spinor));
      
      /* initalize the random sources: one source for each Dirac index beta=src_d --> spin dilution */
      const unsigned int seed_i = measurement_list[id].seed; // has to be same seed
      for (size_t i_phi = 0; i_phi < 2; i_phi++) {    // doublet: light or heavy
        for (size_t src_d = 0; src_d < 4; src_d++) {  // spin dilution index
          for (size_t i_f = 0; i_f < 2; i_f++) {      // flavor index of the doublet
            // light doublet
            eo_source_spinor_field_spin_diluted_oet_ts(arr_eo_spinor[0][i_phi][src_d][0][i_f],
                                                       arr_eo_spinor[0][i_phi][src_d][1][i_f], t0,
                                                       src_d, sample, traj, seed_i);
          }
        }
      }

      // heavy doublet: for each even-odd block, I multiply by (1 + i*tau_2)/sqrt(2)
      // the basis for the inversion should be the same as for the light
      // --> I will multiply later by the inverse, namely at the end of the inversion
      for (size_t src_d = 0; src_d < 4; src_d++) {  // spin dilution index
        for (size_t i_eo = 0; i_eo < 2; i_eo++) {
            // (c,s) --> [(1+i*tau_2)/sqrt(2)] * (c,s) , stored at temporarely in the propagator
            // spinors (used as dummy spinors)
            mul_one_pm_itau2_and_div_by_sqrt2(arr_eo_spinor[1][1][src_d][i_eo][0],
                                              arr_eo_spinor[1][1][src_d][i_eo][0 + 1],
                                              arr_eo_spinor[0][1][src_d][i_eo][0],
                                              arr_eo_spinor[0][1][src_d][i_eo][0 + 1], +1.0,
                                              VOLUME / 2);
          for (size_t i_f = 0; i_f < 2; i_f++) {
            // assigning the result to the first components (the sources).
            // The propagators will be overwritten with the inversion
            assign(arr_eo_spinor[0][1][src_d][i_eo][i_f],
                   arr_eo_spinor[1][1][src_d][i_eo][i_f], VOLUME / 2);
          }
        }
      }

      /*
      assign the sources and propagator pointes for the operators
      these need to be know when calling the inverter() method
      */
      // psi = (psi[(s,p)][phi][beta][eo][f][x])[alpha][c] // last 3 indices are proper of the spinor struct

      for (size_t src_d = 0; src_d < 4; src_d++) {  // spin dilution index
          /* light doublet */

          // up

          // sorces
          optr1->sr0 = arr_eo_spinor[0][0][src_d][0][0];
          optr1->sr1 = arr_eo_spinor[0][0][src_d][1][0];
          // propagators
          optr1->prop0 = arr_eo_spinor[1][0][src_d][0][0];
          optr1->prop1 = arr_eo_spinor[1][0][src_d][1][0];

          optr1->inverter(i1, 0, 0); // inversion for the up flavor

          // inversion of the light doublet only inverts the up block (the operator is diagonal in flavor)
          // down components in flavor will be empty
          optr1->DownProp = 1;

          // sources
          optr1->sr0 = arr_eo_spinor[0][0][src_d][0][1];
          optr1->sr1 = arr_eo_spinor[0][0][src_d][1][1];
          //propagator
          optr1->prop0 = arr_eo_spinor[1][0][src_d][0][1];
          optr1->prop1 = arr_eo_spinor[1][0][src_d][1][1];

          optr1->inverter(i1, 0, 0); //inversion for the down

          optr1->DownProp = 0; //restoring to default

          /* heavy doublet */
          
          optr1->sr0 = arr_eo_spinor[0][1][src_d][0][0];
          optr1->sr1 = arr_eo_spinor[0][1][src_d][1][0];
          optr1->sr2 = arr_eo_spinor[0][1][src_d][0][1];
          optr1->sr3 = arr_eo_spinor[0][1][src_d][1][1];

          optr1->prop0 = arr_eo_spinor[1][1][src_d][0][0];
          optr1->prop1 = arr_eo_spinor[1][1][src_d][1][0];
          optr1->prop2 = arr_eo_spinor[1][1][src_d][0][1];
          optr1->prop3 = arr_eo_spinor[1][1][src_d][1][1];

          optr2->inverter(i2, 0, 0); // inversion for both flavor components
      }


        // conclude the change of basis for the heavy doublet
      for (size_t src_d = 0; src_d < 4; src_d++) {  // spin dilution index
        for (size_t i_eo = 0; i_eo < 2; i_eo++) {
            // (c,s) --> [(1-i*tau_2)/sqrt(2)] * (c,s) , stored at temporarely in the propagator
            // spinors (used as dummy spinors)
            mul_one_pm_itau2_and_div_by_sqrt2(arr_eo_spinor[1][1][src_d][i_eo][0],
                                              arr_eo_spinor[1][1][src_d][i_eo][1],
                                              arr_eo_spinor[0][1][src_d][i_eo][0],
                                              arr_eo_spinor[0][1][src_d][i_eo][1], -1.0,
                                              VOLUME / 2);
          for (size_t i_f = 0; i_f < 2; i_f++) {
            // assigning the result to the first components (the sources).
            // The propagators will be overwritten with the inversion
            assign(arr_eo_spinor[0][1][src_d][i_eo][i_f],
                   arr_eo_spinor[1][1][src_d][i_eo][i_f], VOLUME / 2);
          }
        }
      }

      // (no even-odd): source+propagator, doublet, spin dilution index, even-odd, flavor, position (+ Dirac, color)
      // (psi[i_sp][phi][beta][f][x])[alpha][c] // last 3 indices are proper of the spinor struct
      spinor *****arr_spinor = (spinor *****) callocMultiDimensional({2,2,4,2,VOLUME}, 5, sizeof(spinor));

      // now we switch from even-odd representation to standard
      for (size_t i_sp = 0; i_sp < 2; i_sp++) {         // source or sink
        for (size_t i_phi = 0; i_phi < 2; i_phi++) {    // doublet index
          for (size_t src_d = 0; src_d < 4; src_d++) {  // spin dilution index
            for (size_t i_f = 0; i_f < 2; i_f++) {
              convert_eo_to_lexic(arr_spinor[i_sp][i_phi][src_d][i_f],
                                  arr_eo_spinor[i_sp][i_phi][src_d][0][i_f],
                                  arr_eo_spinor[i_sp][i_phi][src_d][1][i_f]);
            }
          }
        }
      }
      free(arr_eo_spinor); // freeing up space: eo spinors not needed anymore

      /*
        Now that I have all the propagators (all in the basis of
        https://arxiv.org/pdf/1005.2042.pdf) I can build the correlators of eq. (20)
      */

      /* now we sum only over local space for every t */
      for (t = 0; t < T; t++) {
        j = g_ipt[t][0][0][0];

        // dummy variables
        res = 0.;
        respa = 0.;
        resp4 = 0.;

        for (i = j; i < j + LX * LY * LZ; i++) {
          // light correlators
          for (size_t beta = 0; beta < 4; beta++) {  // spin dilution
            double sign = 1.0; // the light source sould be multiplied by gamma_5
            if(beta >= 2){
              sign = -1.0;
            }
            res += sign*_spinor_prod_re(arr_spinor[1][0][beta][0][i], arr_spinor[1][0][beta][0][i]);
            _gamma0(phi, arr_spinor[1][0][beta][0][i]);
            respa += sign*_spinor_prod_re(arr_spinor[1][0][beta][0][i], phi);
            _gamma5(phi, phi);
            resp4 += sign*_spinor_prod_im(arr_spinor[1][0][beta][0][i], phi);
          }

          for (size_t i_f = 0; i_f < 2; i_f++) {
            for (size_t j_f = 0; j_f < 2; j_f++) {
              for (size_t g1 = 0; g1 < 2; g1++) {
                for (size_t g2 = 0; g2 < 2; g2++) {
                  for (size_t beta1 = 0; beta1 < 4; beta1++) {
                    for (size_t beta2 = 0; beta2 < 4; beta2++) {
                      double dum1, dum2 = 0.0, 0.0;  // dummy variables

                      // 1st  contribution
                      phi = arr_spinor[1][1][beta_2][i_f][i];
                      if (Gamma_1 == 1) {  // Gamma_1 = 1
                        _gamma5(phi, phi);
                      }
                      _spinor_scalar_prod(dum1, arr_spinor[0][0][beta_1][1 - i_f][i], phi);

                      // 2nd contribution
                      phi = arr_spinor[1][0][beta_1][1 - j_f][i];
                      if (Gamma_2 == 1) {  // Gamma_1 = gamma_5
                        _gamma5(phi, phi);
                      }
                      _spinor_scalar_prod(dum2, arr_spinor[0][1][beta_2][j_f][i], phi);

                      res_hihj_g1g2[i_f][j_f][beta_1][beta_2] += dum1 * dum2;
                    }
                  }
                }
              }
            }
          }
        }

#if defined TM_USE_MPI
        MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
        res = mpi_res;
        MPI_Reduce(&respa, &mpi_respa, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
        respa = mpi_respa;
        MPI_Reduce(&resp4, &mpi_resp4, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
        resp4 = mpi_resp4;
        sCpp[t] = +res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. / optr1->kappa /
                  optr1->kappa;
        sCpa[t] = -respa / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                  optr1->kappa / optr1->kappa;
        sCp4[t] = +resp4 / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                  optr1->kappa / optr1->kappa;
#else
        Cpp[t] = +res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. / optr1->kappa /
                 optr1->kappa;
        Cpa[t] = -respa / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                 optr1->kappa / optr1->kappa;
        Cp4[t] = +resp4 / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                 optr1->kappa / optr1->kappa;
#endif

        for (size_t hi = 0; hi < 2; hi++) {
          for (size_t hj = 0; hj < 2; hj++) {
            for (size_t g1 = 0; g1 < 2; g1++) {
              for (size_t g2 = 0; g2 < 2; g2++) {
#if defined TM_USE_MPI
                MPI_Reduce(&res_hihj_g1g2[hi][hj][g1][g2], &mpi_res_hihj_g1g2[hi][hj][g1][g2], 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
                res_hihj_g1g2[hi][hj][g1][g2] = mpi_res_hihj_g1g2[hi][hj][g1][g2];
                C_hihj_g1g2[hi][hj][g1][g2][t] = -eta_Gamma[g1] * res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ);
#else
                C_hihj_g1g2[hi][hj][g1][g2][t] = -eta_Gamma[g1] * res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ);
#endif
              }
            }
          }
        }

      } //end loop over "t"

#ifdef TM_USE_MPI
      /* some gymnastics needed in case of parallelisation */
      if (g_mpi_time_rank == 0) {
        // light correlators
        MPI_Gather(sCpp, T, MPI_DOUBLE, Cpp, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCpa, T, MPI_DOUBLE, Cpa, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCp4, T, MPI_DOUBLE, Cp4, T, MPI_DOUBLE, 0, g_mpi_SV_slices);

        // heavy mesons
        for (size_t hi = 0; hi < 2; hi++) {
          for (size_t hj = 0; hj < 2; hj++) {
            for (size_t g1 = 0; g1 < 2; g1++) {
              for (size_t g2 = 0; g2 < 2; g2++) {
                MPI_Gather(sC_hihj_g1g2[hi][hj][g1][g2], T, MPI_DOUBLE, C_hihj_g1g2[hi][hj][g1][g2], T, MPI_DOUBLE, 0, g_mpi_SV_slices);
              }
            }
          }
        }
      }
#endif

      /* and write everything into a file */
      if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
        ofs = fopen(filename, "w");
        fprintf(ofs, "1  1  0  %e  %e\n", Cpp[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "1  1  %d  %e  ", t, Cpp[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%e\n", Cpp[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "1  1  %d  %e  %e\n", t, Cpp[tt], 0.);

        fprintf(ofs, "2  1  0  %e  %e\n", Cpa[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "2  1  %d  %e  ", t, Cpa[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%e\n", Cpa[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "2  1  %d  %e  %e\n", t, Cpa[tt], 0.);

        fprintf(ofs, "6  1  0  %e  %e\n", Cp4[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "6  1  %d  %e  ", t, Cp4[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%e\n", Cp4[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "6  1  %d  %e  %e\n", t, Cp4[tt], 0.);


          // heavy meson correlators
          for (size_t hi = 0; hi < 2; hi++) {
            for (size_t hj = 0; hj < 2; hj++) {
              for (size_t g1 = 0; g1 < 2; g1++) {
                for (size_t g2 = 0; g2 < 2; g2++) {
                   fprintf(ofs, "%d  %d  %d  %d 0  %e  %e\n", hi, hj, g1, g2, C_hihj_g1g2[hi][hj][g1][g1][t0], 0.);
                    for (t = 1; t < g_nproc_t * T / 2; t++) {
                      tt = (t0 + t) % (g_nproc_t * T);
                      fprintf(ofs, "%d  %d  %d  %d  %d  %e  ", hi, hj, g1, g2, t, C_hihj_g1g2[hi][hj][g1][g2][tt]);
                      tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
                      fprintf(ofs, "%e\n", C_hihj_g1g2[hi][hj][g1][g2][tt]);
                    }
                    tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
                    fprintf(ofs, "%d  %d  %d  %d  %d  %d  %e  %e\n", hi, hj, g1, g2, t, C_hihj_g1g2[hi][hj][g1][g2][tt], 0.0);
                }
              }
            }
          }

        fclose(ofs);
      }
#ifdef TM_USE_MPI
      if (g_mpi_time_rank == 0) {
        free(Cpp);
        free(Cpa);
        free(Cp4);
        free(C_hihj_g1g2);
      }
      free(sCpp);
      free(sCpa);
      free(sCp4);
      free(sC_hihj_g1g2);
#else
      free(Cpp);
      free(Cpa);
      free(Cp4);
      free(C_hihj_g1g2);
#endif
    }  // for(max_time_slices)
  }    // for(max_samples)

  tm_stopwatch_pop(&g_timers, 0, 1, "");

  return;
}

void correlators_measurement(const int traj, const int id, const int ieo) {
  light_correlators_measurement(traj, id, ieo);

    // ??? maybe add a double check?
  if (measurement_list[id].measure_heavy_mesons == 1) {
    const unsigned int i1 = 0, i2 = 1; 
    heavy_correlators_measurement(traj, id, ieo, i1, i2);
  }
}
