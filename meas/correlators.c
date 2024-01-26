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
        fprintf( ofs, "1  1  0  %e  %e\n", Cpp[t0], 0.0);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "1  1  %d  %e  ", t, Cpp[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cpp[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "1  1  %d  %e  %e\n", t, Cpp[tt], 0.0);

        fprintf( ofs, "2  1  0  %e  %e\n", Cpa[t0], 0.0);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "2  1  %d  %e  ", t, Cpa[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cpa[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "2  1  %d  %e  %e\n", t, Cpa[tt], 0.0);

        fprintf( ofs, "6  1  0  %e  %e\n", Cp4[t0], 0.0);
        for(t = 1; t < g_nproc_t*T/2; t++) {
          tt = (t0+t)%(g_nproc_t*T);
          fprintf( ofs, "6  1  %d  %e  ", t, Cp4[tt]);
          tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
          fprintf( ofs, "%e\n", Cp4[tt]);
        }
        tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
        fprintf( ofs, "6  1  %d  %e  %e\n", t, Cp4[tt], 0.0);
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


void spinor_dirac_array(su3_vector* phi, spinor psi){
	phi[0] = psi.s0;
	phi[1] = psi.s1;
	phi[2] = psi.s2;
	phi[3] = psi.s3;
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
	printf("Called: %s\n", __func__);

	spinor phi;           // dummy spinor
	// dummy spinors for the heavy doublet inversion
	spinor* phi1 = calloc(VOLUME/2, sizeof(spinor));
	spinor* phi2 = calloc(VOLUME/2, sizeof(spinor));
	int i, j, t, tt, t0;  // dummy indices
	float tmp;            // dummy variable
	
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
	
	/* heavy-light correlators variables */
	// sign change bilinear^\dagger, with Gamma = 1,gamma_5
	double eta_Gamma[2] = {1.0, -1.0};

	// even-odd spinor fields for the light and heavy doublet correlators
	// INDICES: source+propagator, doublet, spin dilution index, flavor projection, even-odd,
	// flavor, position
	// (+ Dirac, color) psi = (psi[(s,p)][db][beta][F][eo][f][x])[alpha][c]
	// last 2 indices come from spinor struct
	// Note: propagator in the sense that it is D^{-1}*source after the inversion
	spinor* arr_eo_spinor[2][2][4][2][2][2];

	// (no even-odd): source+propagator, doublet, spin dilution index, flavor proj, flavor
	// index, position (+ Dirac, color) (psi[i_sp][d][beta][F][f][x])[alpha][c]
	spinor* arr_spinor[2][2][4][2][2];

	// allocating memory
	for (size_t i_sp = 0; i_sp < 2; i_sp++) { // source or propagator
		for (size_t db = 0; db < 2; db++) { // doublet: light or heavy
			for (size_t beta = 0; beta < 4; beta++) { // spin dilution index
				for (size_t F = 0; F < 2; F++) { // flavor dilution index
					for (size_t i_eo = 0; i_eo < 2; i_eo++) { // even-odd index
						for (size_t i_f = 0; i_f < 2; i_f++){// flavor index of the doublet
							arr_eo_spinor[i_sp][db][beta][F][i_eo][i_f] = (spinor*) calloc(VOLUME/2, sizeof(spinor));
							zero_spinor_field(arr_eo_spinor[i_sp][db][beta][F][i_eo][i_f], VOLUME/2);
							
							if (i_eo == 0){ // (trick) doing it only once -> no eo index
								arr_spinor[i_sp][db][beta][F][i_f] = (spinor*) calloc(VOLUME, sizeof(spinor));
								zero_spinor_field(arr_spinor[i_sp][db][beta][F][i_f], VOLUME);
							}
						}
					}
				}
			}
		}
	}
	
	if (g_proc_id == 0){
		printf("Checkpoint 1\n");
	}
	
	if (g_proc_id  == 0){
		printf("Checkpoint 2\n");
	}
	
	double res_hihj_g1g2[2][2][2][2];
	double mpi_res_hihj_g1g2[2][2][2][2];
	
	if (g_proc_id  == 0){
		printf("Checkpoint 3\n");
	}
	
	FILE *ofs;                                // output file stream
	char filename[TM_OMEAS_FILENAME_LENGTH];  // file path
	
	/* building the stochastic propagator spinors needed for the correlators */
	
	// MPI_Barrier(MPI_COMM_WORLD);
	init_operators();  // initialize operators (if not already done)

	if (g_proc_id  == 0){
		printf("Checkpoint 4\n");
	}

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
	operator* optr1 = & operator_list[i1];  // light doublet
	operator* optr2 = & operator_list[i2];  // heavy c,s doublet
	printf("Initialized the operators\n");

	bool b1 = (optr1->type != TMWILSON && optr1->type != WILSON && optr1->type != CLOVER);
	bool b2 = (optr2->type != DBTMWILSON && optr2->type != DBCLOVER);
	if (b1 || b2) {
		if (g_proc_id == 0) {
			if (b1) {
				fprintf(stderr,
								"Warning! optr1 should be one of the following: TMWILSON, WILSON and CLOVER\n");
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
	int max_time_slices = measurement_list[id].all_time_slices ? measurement_list[id].max_source_slice : 1;
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

			// the number of independent correlators is 16 = (2*2)_{h_1 h_2} * (2*2)_{Gamma_1 Gamma_2}
			// hi: c,s
			// Gamma = 1, gamma_5
			double* C_hihj_g1g2[2][2][2][2];
			double* sC_hihj_g1g2[2][2][2][2];
			
			// allocating memory
			for (int h1=0; h1<2; h1++){
				for (int h2=0; h2<2; h2++){
					for (int g1=0; g1<2; g1++){
						for (int g2=0; g2<2; g2++){
#ifdef TM_USE_MPI
							mpi_res_hihj_g1g2[h1][h2][g1][g2] = 0.0;
							res_hihj_g1g2[h1][h2][g1][g2] = 0.0;
							sC_hihj_g1g2[h1][h2][g1][g2] = (double*) calloc(T, sizeof(double));
							if (g_mpi_time_rank == 0) {
								C_hihj_g1g2[h1][h2][g1][g2] = (double*) calloc(g_nproc_t*T, sizeof(double));
							}
#else
							res_hihj_g1g2[h1][h2][g1][g2] = 0.0;
							C_hihj_g1g2[h1][h2][g1][g2] = (double*) calloc(T, sizeof(double));
#endif
						}
					}
				}
			}
			
			/* init random sources: one for each Dirac index beta --> spin dilution */
			const unsigned int seed_i = measurement_list[id].seed;  // has to be same seed
			for (size_t db = 0; db < 2; db++) {                     // doublet: light or heavy
				for (size_t beta = 0; beta < 4; beta++) {             // spin dilution index
					for (size_t F = 0; F < 2; F++) {                    // flavor dilution index
						for (size_t i_f = 0; i_f < 2; i_f++) {            // flavor index of the doublet
							eo_source_spinor_field_spin_diluted_oet_ts(
								arr_eo_spinor[0][db][beta][F][0][i_f],
								arr_eo_spinor[0][db][beta][F][1][i_f], t0,
								beta, sample, traj, seed_i);
						}
					}
				}
			}
			zero_spinor_field(phi1, VOLUME/2);
			zero_spinor_field(phi2, VOLUME/2);
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 5\n");
			}
			
			// heavy doublet: for each even-odd block, I multiply by (1 + i*tau_2)/sqrt(2)
			// the basis for the inversion should be the same as for the light
			// --> I will multiply later by the inverse, namely at the end of the inversion
			for (size_t beta = 0; beta < 4; beta++) {      // spin dilution index
				for (size_t F = 0; F < 2; F++) {             // flavor projection index
					for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd
						// (c,s) --> [(1-i*tau_2)/sqrt(2)] * (c,s)
						// stored temporarely in the propagator spinors (used as dummy)
						if (F==0){
							mul_one_pm_itau2_and_div_by_sqrt2(
							arr_eo_spinor[1][1][beta][F][i_eo][0], arr_eo_spinor[1][1][beta][F][i_eo][1],
							arr_eo_spinor[0][1][beta][F][i_eo][0], phi1, +1.0,
							VOLUME / 2);
						}
						if (F==1){
							mul_one_pm_itau2_and_div_by_sqrt2(
							arr_eo_spinor[1][1][beta][F][i_eo][0], arr_eo_spinor[1][1][beta][F][i_eo][1],
							phi2, arr_eo_spinor[0][1][beta][F][i_eo][1], +1.0,
							VOLUME / 2);
						}
						
						for (size_t i_f = 0; i_f < 2; i_f++) {
							// assigning the result to the first components (the sources).
							// The propagators will be overwritten with the inversion
							assign(arr_eo_spinor[0][1][beta][F][i_eo][i_f],
										 arr_eo_spinor[1][1][beta][F][i_eo][i_f], VOLUME / 2);
							zero_spinor_field(arr_eo_spinor[1][1][beta][F][i_eo][i_f], VOLUME / 2);
						}
					}
				}
			}

			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 6\n");
			}

			/*
				assign the sources and propagator pointes for the operators
				these need to be know when calling the inverter() method */
			for (size_t beta = 0; beta < 4; beta++) {  // spin dilution index
				/* light doublet */
				
				// up
				
				// sources
				optr1->sr0 = arr_eo_spinor[0][0][beta][0][0][0];
				optr1->sr1 = arr_eo_spinor[0][0][beta][0][1][0];
				
				MPI_Barrier(MPI_COMM_WORLD);
				if (g_proc_id == 0){
					printf("Checkpoint 7\n");
				}
				
				// propagators
				optr1->prop0 = arr_eo_spinor[1][0][beta][0][0][0];
				optr1->prop1 = arr_eo_spinor[1][0][beta][0][1][0];
				
				printf("inverting the light doublet\n");
				optr1->inverter(i1, 0, 0);  // inversion for the up flavor
				
				// PLEASE KEEP THESE LINES COMMENTED, MAY BE USEFUL IN THE FUTURE
				// // inversion of the light doublet only inverts the up block (the operator is
				// // diagonal in flavor) down components in flavor will be empty
				// optr1->DownProp = 1;
				
				// // sources
				// optr1->sr0 = arr_eo_spinor[0][0][beta][0][1];
				// optr1->sr1 = arr_eo_spinor[0][0][beta][1][1];
				// // propagator
				// optr1->prop0 = arr_eo_spinor[1][0][beta][0][1];
				// optr1->prop1 = arr_eo_spinor[1][0][beta][1][1];
				
				// optr1->inverter(i1, 0, 0);  // inversion for the down
				
				// optr1->DownProp = 0;        // restoring to default
				
				/* heavy doublet */
				printf("Inverting the heavy doublet\n");
				for (size_t F = 0; F < 2; F++) {             // flavor projection index
					printf("F=%d\n", F);
					
					optr2->sr0 = arr_eo_spinor[0][1][beta][F][0][0];
					optr2->sr1 = arr_eo_spinor[0][1][beta][F][1][0];
					optr2->sr2 = arr_eo_spinor[0][1][beta][F][0][1];
					optr2->sr3 = arr_eo_spinor[0][1][beta][F][1][1];
					
					optr2->prop0 = arr_eo_spinor[1][1][beta][F][0][0];
					optr2->prop1 = arr_eo_spinor[1][1][beta][F][1][0];
					optr2->prop2 = arr_eo_spinor[1][1][beta][F][0][1];
					optr2->prop3 = arr_eo_spinor[1][1][beta][F][1][1];
					
					SourceInfo.no_flavours = 1; // only for the heavy
					optr2->inverter(i2, 0, 0);  // inversion for both flavor components
					SourceInfo.no_flavours = 0; // reset to previous value
				}
			}
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 8\n");
			}
			
			// conclude the change of basis for the heavy doublet
			for (size_t beta = 0; beta < 4; beta++) {      // spin dilution index
				for (size_t F = 0; F < 2; F++) {             // flavor projection index
					for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd
						// (c,s) --> [(1-i*tau_2)/sqrt(2)] * (c,s)
						// stored temporarely in phi1, phi2
						mul_one_pm_itau2_and_div_by_sqrt2(phi1, phi2,
							arr_eo_spinor[1][1][beta][F][i_eo][0], arr_eo_spinor[1][1][beta][F][i_eo][1], -1.0,
							VOLUME / 2);

						assign(arr_eo_spinor[1][1][beta][F][i_eo][0], phi1, VOLUME / 2);
						assign(arr_eo_spinor[1][1][beta][F][i_eo][1], phi2, VOLUME / 2);
					}
				}
			}
			// reset to zero the dummy spinors
			zero_spinor_field(phi1, VOLUME/2);
			zero_spinor_field(phi2, VOLUME/2);
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 9\n");
			}
			
			// now we switch from even-odd representation to standard
			for (size_t i_sp = 0; i_sp < 2; i_sp++) {      // source or sink
				for (size_t db = 0; db < 2; db++) {          // doublet: light of heavy
					for (size_t beta = 0; beta < 4; beta++) {  // spin dilution index
						for (size_t F = 0; F < 2; F++) {         // flavor projection index
							for (size_t i_f = 0; i_f < 2; i_f++) { //  flavor index
								convert_eo_to_lexic(arr_spinor[i_sp][db][beta][F][i_f],
									arr_eo_spinor[i_sp][db][beta][F][0][i_f],
									arr_eo_spinor[i_sp][db][beta][F][1][i_f]);
							}
						}
					}
				}
			}
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 10\n");
			}
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 11.0\n");
			}
			
			/*
				Now that I have all the propagators (all in the basis of
				https://arxiv.org/pdf/1005.2042.pdf) I can build the correlators of eq. (20)
			*/
			const int f0 = 0;  // flavor index of the up

			/* now we sum only over local space for every t */
			const int j_ts = t0-g_proc_coords[0]*T; // checkerboard index of the time at the source

			for (t = 0; t < T; t++) {
				j = g_ipt[t][0][0][0]; // source and sink separated by "t"

				// dummy variables
				res = 0.0;
				respa = 0.0;
				resp4 = 0.0;

				for (i = j; i < j + LX * LY * LZ; i++) {

					// light correlators
					for (size_t beta = 0; beta < 4; beta++) {  // spin dilution
						spinor psi_u = arr_spinor[1][0][beta][f0][f0][i];
						
						res += _spinor_prod_re(psi_u, psi_u);
						_gamma0(phi, psi_u);
						respa += _spinor_prod_re(psi_u, phi);
						_gamma5(phi, phi);
						resp4 += _spinor_prod_im(psi_u, phi);
					}
					
					// heavy correlators						
					for (size_t hi = 0; hi < 2; hi++) {
						for (size_t hj = 0; hj < 2; hj++) {
							for (size_t g1 = 0; g1 < 2; g1++) {
								for (size_t g2 = 0; g2 < 2; g2++) {
									double complex dum_tot = 0.0;
									for (int alpha_1 = 0; alpha_1 < 4; alpha_1++){
										for (int alpha_2 = 0; alpha_2 < 4; alpha_2++){
											// S_\ell
											spinor psi_h = arr_spinor[1][1][alpha_1][hj][hi][i];
											if (g1 == 0){ // Gamma_1 = Id
												_gamma5(psi_h, psi_h);
											}
											su3_vector psi_h_su3[4];
											spinor_dirac_array(&psi_h_su3[0], psi_h);
											spinor psi_u_star = arr_spinor[1][0][alpha_2][f0][f0][i];
											if (g2 == 0){ // Gamma_2 = Id. NOTE: works because Gamma_2=Gamma_2* for Gamma_2=1,gamma_5
												_gamma5(psi_u_star, psi_u_star);
											}
											su3_vector psi_u_star_su3[4];
											spinor_dirac_array(&psi_u_star_su3[0], psi_u_star);
											complex double dum_12 = 0.0;
											_colorvec_scalar_prod(dum_12, psi_u_star_su3[alpha_1], psi_h_su3[alpha_2]);
											dum_tot =+ dum_12;
										}
									}
									res_hihj_g1g2[hi][hj][g1][g2] = creal(dum_tot); // correlators are real
								}
							}
						}
						
					}
				}
				
				const double vol_fact = (g_nproc_x * LX) * (g_nproc_y * LY) * (g_nproc_z * LZ);
#if defined TM_USE_MPI
				MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
				res = mpi_res;
				MPI_Reduce(&respa, &mpi_respa, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
				respa = mpi_respa;
				MPI_Reduce(&resp4, &mpi_resp4, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
				resp4 = mpi_resp4;
				
				sCpp[t] = +res / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
				sCpa[t] = -respa / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
				sCp4[t] = +resp4 / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
#else
				Cpp[t] = +res / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
				Cpa[t] = -respa / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
				Cp4[t] = +resp4 / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
#endif
				
				printf("%d: Checkpoint 11.5\n", g_proc_id);
				
				for (size_t hi = 0; hi < 2; hi++) {
					for (size_t hj = 0; hj < 2; hj++) {
						for (size_t g1 = 0; g1 < 2; g1++) {
							for (size_t g2 = 0; g2 < 2; g2++) {
#if defined TM_USE_MPI
								MPI_Reduce(&res_hihj_g1g2[hi][hj][g1][g2], &mpi_res_hihj_g1g2[hi][hj][g1][g2], 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
								res_hihj_g1g2[hi][hj][g1][g2] = mpi_res_hihj_g1g2[hi][hj][g1][g2];
								
								sC_hihj_g1g2[hi][hj][g1][g2][t] = - eta_Gamma[g1] * res_hihj_g1g2[hi][hj][g1][g2] / vol_fact;
#else
								C_hihj_g1g2[hi][hj][g1][g2][t] = - eta_Gamma[g1] * res_hihj_g1g2[hi][hj][g1][g2] / vol_fact;
#endif
							}
						}
					}
				}
				
			} // end 1st loop over "t"

			
			if (g_proc_id == 0){
				printf("Checkpoint 11.6\n");
			}
			
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
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 12\n");
			}
			
			printf("check: %d, %d\n", g_mpi_time_rank, g_proc_coords[0]);
			
			/* and write everything into a file */
			if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
				printf("Checkpoint 12.1 : -%s-\n", filename);
				ofs = fopen(filename, "w");
				printf("Checkpoint 12.2\n");

				
				fprintf(ofs, "1  1  0  %e  %e\n", Cpp[t0], 0.0);
				for (t = 1; t < g_nproc_t * T / 2; t++) {
					tt = (t0 + t) % (g_nproc_t * T);
					fprintf(ofs, "1  1  %d  %e  ", t, Cpp[tt]);
					tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
					fprintf(ofs, "%e\n", Cpp[tt]);
				}
				tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
				fprintf(ofs, "1  1  %d  %e  %e\n", t, Cpp[tt], 0.0);
				
				fprintf(ofs, "2  1  0  %e  %e\n", Cpa[t0], 0.0);
				for (t = 1; t < g_nproc_t * T / 2; t++) {
					tt = (t0 + t) % (g_nproc_t * T);
					fprintf(ofs, "2  1  %d  %e  ", t, Cpa[tt]);
					tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
					fprintf(ofs, "%e\n", Cpa[tt]);
				}
				tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
				fprintf(ofs, "2  1  %d  %e  %e\n", t, Cpa[tt], 0.0);

				fprintf(ofs, "6  1  0  %e  %e\n", Cp4[t0], 0.0);
				for (t = 1; t < g_nproc_t * T / 2; t++) {
					tt = (t0 + t) % (g_nproc_t * T);
					fprintf(ofs, "6  1  %d  %e  ", t, Cp4[tt]);
					tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
					fprintf(ofs, "%e\n", Cp4[tt]);
				}
				tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
				fprintf(ofs, "6  1  %d  %e  %e\n", t, Cp4[tt], 0.0);
			 
				printf("Checkpoint 12.3\n");
				// heavy mesons correlators
				for (size_t hi = 0; hi < 2; hi++) {
					for (size_t hj = 0; hj < 2; hj++) {
						for (size_t g1 = 0; g1 < 2; g1++) {
							for (size_t g2 = 0; g2 < 2; g2++) {
								fprintf(ofs, "%u  %u  %u  %u  0  %e  %e\n", hi, hj, g1, g2, C_hihj_g1g2[hi][hj][g1][g1][t0], 0.0);
								for (t = 1; t < g_nproc_t * T / 2; t++) {
									tt = (t0 + t) % (g_nproc_t * T);
									fprintf(ofs, "%u  %u  %u  %u  %d  %e  ", hi, hj, g1, g2, t, C_hihj_g1g2[hi][hj][g1][g2][tt]);
									tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
									fprintf(ofs, "%e\n", C_hihj_g1g2[hi][hj][g1][g2][tt]);
								}
								tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
								fprintf(ofs, "%u  %u  %u  %u  %d  %e  %e\n", hi, hj, g1, g2, t, C_hihj_g1g2[hi][hj][g1][g2][tt], 0.0);
							}
						}
					}
				}
				fclose(ofs);
			}
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 13\n");
			}
			
			// freeing memory: light correlators
#ifdef TM_USE_MPI
			if (g_mpi_time_rank == 0) {
				free(Cpp);
				free(Cpa);
				free(Cp4);
			}
			free(sCpp);
			free(sCpa);
			free(sCp4);
#else
			free(Cpp);
			free(Cpa);
			free(Cp4);
#endif
			
			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 14\n");
			}

			
			// freeing memory: heavy-light correlators
			for (int h1=0; h1<2; h1++){
				for (int h2=0; h2<2; h2++){
					for (int g1=0; g1<2; g1++){
						for (int g2=0; g2<2; g2++){
#ifdef TM_USE_MPI
							free(sC_hihj_g1g2[h1][h2][g1][g2]);
							if (g_mpi_time_rank == 0) {
								free(C_hihj_g1g2[h1][h2][g1][g2]);
							}
#else
							free(C_hihj_g1g2[h1][h2][g1][g2]);
#endif
						}
					}
				}
			}

			//MPI_Barrier(MPI_COMM_WORLD);
			if (g_proc_id == 0){
				printf("Checkpoint 15\n");
			}
			
		}  // for(max_time_slices)
	}    // for(max_samples)
	
	// freeing up memory
	for (size_t i_sp = 0; i_sp < 2; i_sp++) { // source or propagator
		for (size_t db = 0; db < 2; db++) { // doublet: light or heavy
			for (size_t beta = 0; beta < 4; beta++) { // spin dilution index
				for (size_t F = 0; F < 2; F++) { // flavor dilution index
					for (size_t i_eo = 0; i_eo < 2; i_eo++) { // even-odd index
						for (size_t i_f = 0; i_f < 2; i_f++){// flavor index of the doublet
							free(arr_eo_spinor[i_sp][db][beta][F][i_eo][i_f]);
							if(i_eo == 0){ // only once: no "eo" index for arr_spinor
								free(arr_spinor[i_sp][db][beta][F][i_f]);
							}
						}
					}
				}
			}
		}
	}

	free(phi1);
	free(phi2);
	
	tm_stopwatch_pop(&g_timers, 0, 1, "");
	
	return;
}

void correlators_measurement(const int traj, const int id, const int ieo) {
	
	// ??? maybe add a double check? does i1 --> light and i2 --> heavy?
	if (measurement_list[id].measure_heavy_mesons == 1) {
		const int i1 = 0, i2 = 1;
		heavy_correlators_measurement(traj, id, ieo, i1, i2);
	}

	light_correlators_measurement(traj, id, ieo);
}
