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
#include <tmlqcd_config.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "correlators.h"
#include "geometry_eo.h"
#include "gettime.h"
#include "global.h"
#include "invert_eo.h"
#include "io/params.h"
#include "io/spinor.h"
#include "io/utils.h"
#include "linalg/convert_eo_to_lexic.h"
#include "measurements.h"
#include "operator.h"
#include "operator/tm_operators_nd.h"
#include "ranlxs.h"
#include "solver/solver.h"
#include "source_generation.h"
#include "start.h"
#include "su3spinor.h"

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
  // int i, j, t, tt, t0;
  double *Cpp = NULL, *Cpa = NULL, *Cp4 = NULL;
  double res = 0., respa = 0., resp4 = 0.;

  float tmp;
  operator* optr;
#ifdef TM_USE_MPI
  double mpi_res = 0., mpi_respa = 0., mpi_resp4 = 0.;
  // send buffer for MPI_Gather
  double *sCpp = NULL, *sCpa = NULL, *sCp4 = NULL;
#endif
  FILE *ofs;
  char filename[TM_OMEAS_FILENAME_LENGTH];
  spinor phi;

  init_operators();
  if (no_operators < 1) {
    if (g_proc_id == 0) {
      fprintf(stderr,
              "Warning! no operators defined in input file, cannot perform online correlator "
              "mesurements!\n");
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");
    return;
  }
  if (no_operators > 1 && g_proc_id == 0) {
    fprintf(stderr, "Warning! number of operators defined larger than 1, using only the first!\n");
  }
  optr = &operator_list[0];
  // we don't want to do inversion twice for this purpose here
  optr->DownProp = 0;
  if (optr->type != TMWILSON && optr->type != WILSON && optr->type != CLOVER) {
    if (g_proc_id == 0) {
      fprintf(stderr,
              "Warning! correlator online measurement currently only implemented for TMWILSON, "
              "WILSON and CLOVER\n");
      fprintf(stderr, "Cannot perform correlator online measurement!\n");
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");
    return;
  }

  if (ranlxs_init == 0) {
    rlxs_init(1, 123456);
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
      if (max_samples == 1 && max_time_slices == 1) {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s%06d", "onlinemeas.", traj);
      } else if (max_samples == 1 && max_time_slices > 1) {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s.t%03d.%06d", "onlinemeas", ts, traj);
      } else {
        snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "%s.s%03d.%06d", "onlinemeas", sample, traj);
      }
      /* generate random timeslice */
      int t0 = ts;
      if (!measurement_list[id].all_time_slices) {
        ranlxs(&tmp, 1);
        t0 = (int)(measurement_list[id].max_source_slice * tmp);
      }
#ifdef TM_USE_MPI
      MPI_Bcast(&t0, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
      if (g_debug_level > 1 && g_proc_id == 0) {
        printf("# timeslice set to %d (T=%d) for online measurement\n", t0, g_nproc_t * T);
        printf("# online measurements parameters: kappa = %.12f, mu = %.12f\n", optr->kappa,
               optr->mu / 2. / optr->kappa);
      }
      // atime = gettime();

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
      source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], t0, sample, traj,
                                  measurement_list[id].seed);
      optr->sr0 = g_spinor_field[0];
      optr->sr1 = g_spinor_field[1];
      optr->prop0 = g_spinor_field[2];
      optr->prop1 = g_spinor_field[3];
      // the propagator fields are the initial guess of the inversion (e.g. CG): start from zero,
      // not from whatever the HMC or a previous measurement left in these work fields
      zero_spinor_field(g_spinor_field[2], VOLUME / 2);
      zero_spinor_field(g_spinor_field[3], VOLUME / 2);

      // op_id = 0, index_start = 0, write_prop = 0
      optr->inverter(0, 0, 0);

      /* now we bring it to normal format */
      /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
      convert_eo_to_lexic(g_spinor_field[DUM_MATRIX], g_spinor_field[2], g_spinor_field[3]);

      /* now we sum only over local space for every t */
      for (int t = 0; t < T; t++) {
        int j = g_ipt[t][0][0][0];
        res = 0.;
        respa = 0.;
        resp4 = 0.;
        for (int i = j; i < j + LX * LY * LZ; i++) {
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
        sCpp[t] = +res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. / optr->kappa /
                  optr->kappa;
        sCpa[t] = -respa / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                  optr->kappa / optr->kappa;
        sCp4[t] = +resp4 / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                  optr->kappa / optr->kappa;
#else
        Cpp[t] = +res / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. / optr->kappa /
                 optr->kappa;
        Cpa[t] = -respa / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                 optr->kappa / optr->kappa;
        Cp4[t] = +resp4 / (g_nproc_x * LX) / (g_nproc_y * LY) / (g_nproc_z * LZ) / 2. /
                 optr->kappa / optr->kappa;
#endif
      }

#ifdef TM_USE_MPI
      /* some gymnastics needed in case of parallelisation */
      if (g_mpi_time_rank == 0) {
        MPI_Gather(sCpp, T, MPI_DOUBLE, Cpp, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCpa, T, MPI_DOUBLE, Cpa, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
        MPI_Gather(sCp4, T, MPI_DOUBLE, Cp4, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
      }
#endif

      /* and write everything into a file */
      if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
        int tt, t;
        ofs = fopen(filename, "w");
        fprintf(ofs, "1  1  0  %.16e  %.16e\n", Cpp[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "1  1  %d  %.16e  ", t, Cpp[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%.16e\n", Cpp[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "1  1  %d  %.16e  %.16e\n", t, Cpp[tt], 0.);

        fprintf(ofs, "2  1  0  %.16e  %.16e\n", Cpa[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "2  1  %d  %.16e  ", t, Cpa[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%.16e\n", Cpa[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "2  1  %d  %.16e  %.16e\n", t, Cpa[tt], 0.);

        fprintf(ofs, "6  1  0  %.16e  %.16e\n", Cp4[t0], 0.);
        for (t = 1; t < g_nproc_t * T / 2; t++) {
          tt = (t0 + t) % (g_nproc_t * T);
          fprintf(ofs, "6  1  %d  %.16e  ", t, Cp4[tt]);
          tt = (t0 + g_nproc_t * T - t) % (g_nproc_t * T);
          fprintf(ofs, "%.16e\n", Cp4[tt]);
        }
        tt = (t0 + g_nproc_t * T / 2) % (g_nproc_t * T);
        fprintf(ofs, "6  1  %d  %.16e  %.16e\n", t, Cp4[tt], 0.);
        fclose(ofs);
      }
#ifdef TM_USE_MPI
      if (g_mpi_time_rank == 0) {
        free(Cpp);
        free(Cpa);
        free(Cp4);
      }
      free(sCpp);
      free(sCpa);
      free(sCp4);
      MPI_Barrier(MPI_COMM_WORLD);
#else
      free(Cpp);
      free(Cpa);
      free(Cp4);
#endif
    }  // for(max_time_slices)
  }  // for(max_samples)
  // etime = gettime();
  // if(g_proc_id == 0 && g_debug_level > 0) {
  //   printf("ONLINE: measurement done int t/s = %1.4e\n", etime - atime);
  // }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}

// write one even-odd spinor field as a record of a LIME file (collective: all ranks must call it)
// append = 0 --> new file, 1 --> append the record to the existing file
// NOTE: the etmc-source/propagator-format records continue a LIME message (MB=0) and cannot be
// the first record of a writer, so only a self-contained propagator-type record is written
// when the file is created. Each field can be read back with read_spinor(even, odd, filename, n).
static void write_eo_spinor_field(char *filename, spinor *even, spinor *odd, const int append) {
  WRITER *writer = NULL;
  construct_writer(&writer, filename, append);
  if (!append) {
    write_propagator_type(writer, 0);
  }
  write_spinor(writer, &even, &odd, 1, 64);
  destruct_writer(writer);
}

// correlators of the heavy-meson measurement per time slice: 3 light (Cpp, Cpa, Cp4) at 0..2 and
// 16 heavy at hm_heavy_index(hi, hj, g1, g2): 16 = (2*2)_{h_1 h_2} * (2*2)_{Gamma_1 Gamma_2},
// hi: c,s, Gamma = 1, gamma_5
#define HM_LIGHT 3
#define HM_N_CORR (HM_LIGHT + 16)
static inline size_t hm_heavy_index(const size_t hi, const size_t hj, const size_t g1,
                                    const size_t g2) {
  return HM_LIGHT + ((hi * 2 + hj) * 2 + g1) * 2 + g2;
}

// contractions of the heavy-meson measurement for one spin dilution index beta: adds the sums
// over the local spatial sites of every local time slice t to corr_thr[thread][t][k]
// (k: see above), i.e. per-thread partial sums. The signs and normalisations are applied after
// the sums over beta, threads and ranks.
// prop_eo[db][F][eo][f]: propagators of this beta (heavy doublet in the basis of
// https://arxiv.org/pdf/1005.2042.pdf), see heavy_correlators_measurement
static void heavy_mesons_contractions(double *const corr_thr, spinor *prop_eo[2][2][2][2],
                                      const size_t beta) {
  const double g5_diag[4] = {+1.0, +1.0, -1.0, -1.0};  // diagonal of gamma_5 (chiral basis)
#ifdef TM_USE_OMP
#pragma omp parallel
#endif
  {
#ifdef TM_USE_OMP
    double *const acc = corr_thr + (size_t)omp_get_thread_num() * T * HM_N_CORR;
#else
    double *const acc = corr_thr;
#endif
    for (int t = 0; t < T; t++) {
      const int j = g_ipt[t][0][0][0];  // source and sink separated by "t"
      double sum[HM_N_CORR] = {0.0};
      // static schedule: the same partition of the sites for a given number of threads
#ifdef TM_USE_OMP
#pragma omp for schedule(static) nowait
#endif
      for (int i = j; i < j + LX * LY * LZ; i++) {
        // the propagators are read directly from the even-odd fields: parity and index of the
        // lexicographic site i (as in convert_eo_to_lexic)
        const size_t eo = (g_lexic2eo[i] < (VOLUME + RAND) / 2) ? 0 : 1;
        const int ix = g_lexic2eosub[i];
        const spinor *const psi_u = &prop_eo[0][0][eo][0][ix];  // light: up flavor (F = f = 0)

        // light correlators
        spinor phi;
        sum[0] += _spinor_prod_re(*psi_u, *psi_u);
        _gamma0(phi, *psi_u);
        sum[1] += _spinor_prod_re(*psi_u, phi);
        _gamma5(phi, phi);
        sum[2] += _spinor_prod_im(*psi_u, phi);

        // heavy correlators:
        // g_{s2,s1} [S^h]_{s1,be} g_{be,be} [s^u*]_{s2,be}
        // the matrix g can be 1 or gamma_5
        for (size_t hi = 0; hi < 2; hi++) {
          for (size_t hj = 0; hj < 2; hj++) {
            // the scalar product <psi_u, g psi_h> depends only on g1 (g1 = 0: gamma_5 psi_h,
            // g1 = 1: psi_h), g2 only multiplies by g5_diag[beta]
            const spinor *const psi_h = &prop_eo[1][hj][eo][hi][ix];
            spinor g5_psi_h;
            _gamma5(g5_psi_h, *psi_h);
            double complex prod[2];
            _spinor_scalar_prod(prod[0], *psi_u, g5_psi_h);
            _spinor_scalar_prod(prod[1], *psi_u, *psi_h);
            for (size_t g1 = 0; g1 < 2; g1++) {
              for (size_t g2 = 0; g2 < 2; g2++) {
                // mixed g1 != g2: e_g1*e_g2 = i and Re(i*z) = -Im(z) (sign applied later)
                const double v = (g1 == g2) ? creal(prod[g1]) : cimag(prod[g1]);
                sum[hm_heavy_index(hi, hj, g1, g2)] += (g2 == 0) ? g5_diag[beta] * v : v;
              }
            }
          }
        }
      }
      for (size_t k = 0; k < HM_N_CORR; k++) {
        acc[t * HM_N_CORR + k] += sum[k];
      }
    }
  }  // omp parallel
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

  // dummy spinors for the heavy doublet inversion (allocated after the checks below)
  spinor *phi1 = NULL, *phi2 = NULL;
  int t, tt, t0;  // dummy indices
  float tmp;      // dummy variable

  // even-odd spinor fields for the light and heavy doublet correlators
  // INDICES: doublet, flavor projection, even-odd, flavor, position
  // (+ Dirac, color from the spinor struct), e.g. (prop_eo[db][F][eo][f][x])[alpha][c]
  // The correlators are sums over the spin dilution index beta: sources, inversions and
  // contractions are done for one beta at a time, so one set of sources and propagators
  // (D^{-1}*source) is reused for all beta.
  // Unused combinations are not allocated (NULL), see below.
  spinor *src_eo[2][2][2][2];
  spinor *prop_eo[2][2][2][2];

  FILE *ofs;                                // output file stream
  char filename[TM_OMEAS_FILENAME_LENGTH];  // file path

  /* building the stochastic propagator spinors needed for the correlators */

  // MPI_Barrier(MPI_COMM_WORLD);
  init_operators();  // initialize operators (if not already done)

  if (i1 >= no_operators || i2 >= no_operators) {
    if (g_proc_id == 0) {
      // we don't want the simulation to stop, we can do the measurements offline afterwards
      fprintf(stderr,
              "Warning! heavy-meson online measurement needs operators %d (light) and %d (heavy "
              "doublet), but only %d operator(s) defined in input file, cannot perform heavy-meson "
              "online measurement!\n",
              i1, i2, no_operators);
    }
    tm_stopwatch_pop(&g_timers, 0, 0, "");  // timer for profiling
    return;
  }

  // selecting the operators i1 and i2 from the list of operators initialized before
  operator* optr1 = & operator_list[i1];  // light doublet
  operator* optr2 = & operator_list[i2];  // heavy c,s doublet
  // as in light_correlators_measurement: with AddDownPropagator = yes, op_invert would
  // additionally invert at -mu and return that propagator instead of the +mu one
  optr1->DownProp = 0;

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

  // allocating memory (only after the checks above, so that their early returns do not leak)
  // Only the fields that are used are allocated (calloc: zero-initialised), the others are NULL:
  // light doublet: only the up flavor is inverted -> F = 0, f = 0
  phi1 = calloc(VOLUME / 2, sizeof(spinor));
  phi2 = calloc(VOLUME / 2, sizeof(spinor));
  for (size_t db = 0; db < 2; db++) {            // doublet: light or heavy
    for (size_t F = 0; F < 2; F++) {             // flavor dilution index
      for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd index
        for (size_t i_f = 0; i_f < 2; i_f++) {   // flavor index of the doublet
          const bool used = (db == 1 || (F == 0 && i_f == 0));
          // these fields are passed to the inverters (sr*, prop*), whose Hopping_Matrix
          // exchanges the boundary into the field itself without halfspinor
          // -> they need the halo, i.e. VOLUMEPLUSRAND/2
          src_eo[db][F][i_eo][i_f] =
              used ? (spinor *)calloc(VOLUMEPLUSRAND / 2, sizeof(spinor)) : NULL;
          prop_eo[db][F][i_eo][i_f] =
              used ? (spinor *)calloc(VOLUMEPLUSRAND / 2, sizeof(spinor)) : NULL;
        }
      }
    }
  }

  // correlators (k: see hm_heavy_index), summed over beta per thread and local time slice
  // (corr_thr[thread][t][k]), then over the threads in thread order (deterministic for a given
  // number of threads) and with MPI over the ranks of each time slice; corr[t][k] holds all
  // global time slices on the rank that writes the output
#ifdef TM_USE_OMP
  const int nthreads = omp_num_threads;
#else
  const int nthreads = 1;
#endif
  double *corr_thr = calloc((size_t)nthreads * T * HM_N_CORR, sizeof(double));
  double *corr = calloc((size_t)g_nproc_t * T * HM_N_CORR, sizeof(double));
#ifdef TM_USE_MPI
  double *corr_loc = calloc(T * HM_N_CORR, sizeof(double));  // sum over the threads
  double *corr_red = calloc(T * HM_N_CORR, sizeof(double));  // ... and over the ranks of t
#else
  double *corr_loc = corr;  // without MPI the sum over the threads is the result
#endif

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
      tm_stopwatch_push(&g_timers, "single_sample_single_time_slice", "");  // timer for profiling
      // setting output filename (suffix also used for the source and propagator files)
      char suffix[32];
      if (max_samples == 1 && max_time_slices == 1) {
        snprintf(suffix, sizeof(suffix), "%06d", traj);
      } else if (max_samples == 1 && max_time_slices > 1) {
        snprintf(suffix, sizeof(suffix), "t%03d.%06d", ts, traj);
      } else {
        snprintf(suffix, sizeof(suffix), "s%03d.%06d", sample, traj);
      }
      snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "heavy_mesons.%s", suffix);

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

      // the sums over beta start from zero for every sample / time slice
      for (size_t k = 0; k < (size_t)nthreads * T * HM_N_CORR; k++) {
        corr_thr[k] = 0.0;
      }

      /* Sources, inversions and contractions, one spin dilution index beta at a time: the
         correlators are sums over beta, so one set of source and propagator fields is reused.
         Light and heavy doublet need the same noise. The source depends only on
         (t0, beta, sample, traj, seed): it is generated for the light doublet and the heavy
         doublet sources are built from it with the flavor dilution F = 0 -> (eta, 0),
         F = 1 -> (0, eta), where the zero components are taken from phi1, phi2. */
      const unsigned int seed_i = measurement_list[id].seed;  // has to be same seed
      char src_filename[TM_OMEAS_FILENAME_LENGTH];
      snprintf(src_filename, TM_OMEAS_FILENAME_LENGTH, "heavy_mesons_source.%s", suffix);
      // propagator files, the records are appended for every beta:
      // light: 4 records, beta = 0..3 (up flavor), normalised with 2*kappa (as in op_invert)
      // heavy: 16 records, loop order beta = 0..3, F = 0..1 (source flavor), f = 0..1
      //        (sink flavor), basis of https://arxiv.org/pdf/1005.2042.pdf (index 0 = c, 1 = s),
      //        raw inverter normalisation (no 2*kappa factor)
      char prop_filename_light[TM_OMEAS_FILENAME_LENGTH];
      char prop_filename_heavy[TM_OMEAS_FILENAME_LENGTH];
      snprintf(prop_filename_light, TM_OMEAS_FILENAME_LENGTH, "heavy_mesons_prop_light.%s", suffix);
      snprintf(prop_filename_heavy, TM_OMEAS_FILENAME_LENGTH, "heavy_mesons_prop_heavy.%s", suffix);
      for (size_t beta = 0; beta < 4; beta++) {  // spin dilution index
        tm_stopwatch_push(&g_timers, "sources_generation", "");
        eo_source_spinor_field_spin_diluted_oet_ts(src_eo[0][0][0][0], src_eo[0][0][1][0], t0,
                                                   beta, sample, traj, seed_i);
        tm_stopwatch_pop(&g_timers, 0, 1, "");

        // writing the sources: 4 records, one for each spin dilution index beta = 0..3.
        // The noise is the same for the light and heavy doublet (before the tau_2 rotation),
        // so only the light sources are written. Source time slice: t0.
        if (measurement_list[id].write_heavy_mesons_sources) {
          tm_stopwatch_push(&g_timers, "writing_sources", "");
          write_eo_spinor_field(src_filename, src_eo[0][0][0][0], src_eo[0][0][1][0], (beta > 0));
          tm_stopwatch_pop(&g_timers, 0, 1, "");
        }

        tm_stopwatch_push(&g_timers, "sources_rotation", "");
        // phi1, phi2 are the zero components here (and temporary fields of the propagator
        // rotation of the previous beta, see below)
        zero_spinor_field(phi1, VOLUME / 2);
        zero_spinor_field(phi2, VOLUME / 2);
        // heavy doublet: for each even-odd block, I multiply by (1 + i*tau_2)/sqrt(2)
        // the basis for the inversion should be the same as for the light
        // --> I will multiply later by the inverse, namely at the end of the inversion
        for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd
          // (c,s) --> [(1-i*tau_2)/sqrt(2)] * (c,s), with (c,s) = (eta, 0) for F = 0 and
          // (0, eta) for F = 1
          mul_one_pm_itau2(src_eo[1][0][i_eo][0], src_eo[1][0][i_eo][1], src_eo[0][0][i_eo][0],
                           phi1, +1.0, VOLUME / 2);
          mul_one_pm_itau2(src_eo[1][1][i_eo][0], src_eo[1][1][i_eo][1], phi2,
                           src_eo[0][0][i_eo][0], +1.0, VOLUME / 2);
          // the propagator fields are the initial guess of the inversions (e.g. CG): start from
          // zero, not from the solution of the previous beta / sample / time slice
          zero_spinor_field(prop_eo[0][0][i_eo][0], VOLUME / 2);
          for (size_t F = 0; F < 2; F++) {
            for (size_t i_f = 0; i_f < 2; i_f++) {
              zero_spinor_field(prop_eo[1][F][i_eo][i_f], VOLUME / 2);
            }
          }
        }
        tm_stopwatch_pop(&g_timers, 0, 1, "");

        tm_stopwatch_push(&g_timers, "sources_inversion", "");
        /*
                assign the sources and propagator pointes for the operators
                these need to be know when calling the inverter() method */
        /* light doublet */

        // up

        // sources
        optr1->sr0 = src_eo[0][0][0][0];
        optr1->sr1 = src_eo[0][0][1][0];

        // propagators
        optr1->prop0 = prop_eo[0][0][0][0];
        optr1->prop1 = prop_eo[0][0][1][0];

        optr1->inverter(i1, 0, 0);  // inversion for the up flavor

        // PLEASE KEEP THESE LINES COMMENTED, MAY BE USEFUL IN THE FUTURE
        // // inversion of the light doublet only inverts the up block (the operator is
        // // diagonal in flavor) down components in flavor will be empty
        // // (the light f = 1 fields are not allocated, see the allocation above)
        // optr1->DownProp = 1;

        // // sources
        // optr1->sr0 = src_eo[0][0][0][1];
        // optr1->sr1 = src_eo[0][0][1][1];
        // // propagator
        // optr1->prop0 = prop_eo[0][0][0][1];
        // optr1->prop1 = prop_eo[0][0][1][1];

        // optr1->inverter(i1, 0, 0);  // inversion for the down

        // optr1->DownProp = 0;        // restoring to default

        /* heavy doublet */
        for (size_t F = 0; F < 2; F++) {  // flavor projection index

          optr2->sr0 = src_eo[1][F][0][0];
          optr2->sr1 = src_eo[1][F][1][0];
          optr2->sr2 = src_eo[1][F][0][1];
          optr2->sr3 = src_eo[1][F][1][1];

          optr2->prop0 = prop_eo[1][F][0][0];
          optr2->prop1 = prop_eo[1][F][1][0];
          optr2->prop2 = prop_eo[1][F][0][1];
          optr2->prop3 = prop_eo[1][F][1][1];

          const int no_flavours = SourceInfo.no_flavours;
          SourceInfo.no_flavours = 1;            // only for the heavy
          optr2->inverter(i2, 0, 0);             // inversion for both flavor components
          SourceInfo.no_flavours = no_flavours;  // reset to previous value
        }
        tm_stopwatch_pop(&g_timers, 0, 1, "");

        tm_stopwatch_push(&g_timers, "propagators_rotations", "");
        // conclude the change of basis for the heavy doublet
        for (size_t F = 0; F < 2; F++) {             // flavor projection index
          for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd
            // (c,s) --> [(1-i*tau_2)/sqrt(2)] * (c,s)
            // stored temporarely in phi1, phi2
            mul_one_pm_itau2(phi1, phi2, prop_eo[1][F][i_eo][0], prop_eo[1][F][i_eo][1], -1.0,
                             VOLUME / 2);

            assign(prop_eo[1][F][i_eo][0], phi1, VOLUME / 2);
            assign(prop_eo[1][F][i_eo][1], phi2, VOLUME / 2);
          }
        }
        tm_stopwatch_pop(&g_timers, 0, 1, "");

        // writing the propagators of this beta (file contents: see above)
        if (measurement_list[id].write_heavy_mesons_props) {
          tm_stopwatch_push(&g_timers, "writing_propagators", "");
          write_eo_spinor_field(prop_filename_light, prop_eo[0][0][0][0], prop_eo[0][0][1][0],
                                (beta > 0));
          for (size_t F = 0; F < 2; F++) {
            for (size_t i_f = 0; i_f < 2; i_f++) {
              write_eo_spinor_field(prop_filename_heavy, prop_eo[1][F][0][i_f],
                                    prop_eo[1][F][1][i_f], (beta > 0 || F > 0 || i_f > 0));
            }
          }
          tm_stopwatch_pop(&g_timers, 0, 1, "");
        }

        /*
                Now that I have the propagators of this beta (all in the basis of
                https://arxiv.org/pdf/1005.2042.pdf) I can add their contribution to the
                correlators of eq. (20)
        */
        tm_stopwatch_push(&g_timers, "contractions", "");
        heavy_mesons_contractions(corr_thr, prop_eo, beta);
        tm_stopwatch_pop(&g_timers, 0, 1, "");
      }  // beta

      tm_stopwatch_push(&g_timers, "contractions", "");
      // sum over the threads (in thread order)
      for (size_t k = 0; k < (size_t)T * HM_N_CORR; k++) {
        double sum = 0.0;
        for (int thr = 0; thr < nthreads; thr++) {
          sum += corr_thr[(size_t)thr * T * HM_N_CORR + k];
        }
        corr_loc[k] = sum;
      }
#ifdef TM_USE_MPI
      // sum over the ranks of each local time slice and gather of all time slices on the rank
      // that writes the output: one reduction and one gather of all T * HM_N_CORR values
      MPI_Reduce(corr_loc, corr_red, T * HM_N_CORR, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
      if (g_mpi_time_rank == 0) {
        MPI_Gather(corr_red, T * HM_N_CORR, MPI_DOUBLE, corr, T * HM_N_CORR, MPI_DOUBLE, 0,
                   g_mpi_SV_slices);
      }
#endif
      tm_stopwatch_pop(&g_timers, 0, 1, "");

      /* normalise and write everything into a file */
      if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
        const int T_glob = g_nproc_t * T;
        const double vol_fact = (g_nproc_x * LX) * (g_nproc_y * LY) * (g_nproc_z * LZ);
        // the operators are defined in the basis of https://arxiv.org/pdf/1005.2042.pdf
        // Oc1 = -  bar_chi^d Id chi^c    hi=0   g1=0
        // Oc5 = -i bar_chi^d g5 chi^c    hi=0   g1=1
        // Os1 = +  bar_chi^d Id chi^s    hi=1   g1=0
        // Os5 = +i bar_chi^d g5 chi^s    hi=1   g1=1
        // each operator has a prefactor s_h*e_g with s_c=-1, s_s=1 and e_g1=1, e_g5=i
        // the dagger of the operator do not produce any extra sing but there is a -1 fermionic
        // loop sign -1 * s_hi * s_hj * e_g1 * e_g2
        // (the signs are +-1, so applying them after the summation is exact)
        const int sign_hihj[2][2] = {{1, -1}, {-1, 1}};  // sign for the heavy doublet
        // mixed g1 != g2 e_g1*e_g2 = i and Re(i*z) = -Im(z) -> sign -1 with cimag
        const int sign_g1g2[2][2] = {{1, -1}, {-1, -1}};  // sign for the gamma matrices
        for (t = 0; t < T_glob; t++) {
          double *const c = corr + (size_t)t * HM_N_CORR;
          c[0] = +c[0] / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
          c[1] = -c[1] / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
          c[2] = +c[2] / vol_fact / 2.0 / optr1->kappa / optr1->kappa;
          for (size_t hi = 0; hi < 2; hi++) {
            for (size_t hj = 0; hj < 2; hj++) {
              for (size_t g1 = 0; g1 < 2; g1++) {
                for (size_t g2 = 0; g2 < 2; g2++) {
                  const size_t k = hm_heavy_index(hi, hj, g1, g2);
                  c[k] = -1 * (sign_hihj[hi][hj] * sign_g1g2[g1][g2] * c[k]) / vol_fact;
                }
              }
            }
          }
        }

        ofs = fopen(filename, "w");

        // light correlators Cpp, Cpa, Cp4
        const int light_id[HM_LIGHT] = {1, 2, 6};
        for (size_t k = 0; k < HM_LIGHT; k++) {
          fprintf(ofs, "%d  1  0  %.16e  %.16e\n", light_id[k], corr[t0 * HM_N_CORR + k], 0.0);
          for (t = 1; t < T_glob / 2; t++) {
            tt = (t0 + t) % T_glob;
            fprintf(ofs, "%d  1  %d  %.16e  ", light_id[k], t, corr[tt * HM_N_CORR + k]);
            tt = (t0 + T_glob - t) % T_glob;
            fprintf(ofs, "%.16e\n", corr[tt * HM_N_CORR + k]);
          }
          tt = (t0 + T_glob / 2) % T_glob;
          fprintf(ofs, "%d  1  %d  %.16e  %.16e\n", light_id[k], t, corr[tt * HM_N_CORR + k], 0.0);
        }

        // heavy mesons correlators
        for (size_t hi = 0; hi < 2; hi++) {
          for (size_t hj = 0; hj < 2; hj++) {
            for (size_t g1 = 0; g1 < 2; g1++) {
              for (size_t g2 = 0; g2 < 2; g2++) {
                const size_t k = hm_heavy_index(hi, hj, g1, g2);
                fprintf(ofs, "%zu  %zu  %zu  %zu  0  %.16e  %.16e\n", hi, hj, g1, g2,
                        corr[t0 * HM_N_CORR + k], 0.0);
                for (t = 1; t < T_glob / 2; t++) {
                  tt = (t0 + t) % T_glob;
                  fprintf(ofs, "%zu  %zu  %zu  %zu  %d  %.16e  ", hi, hj, g1, g2, t,
                          corr[tt * HM_N_CORR + k]);
                  tt = (t0 + T_glob - t) % T_glob;
                  fprintf(ofs, "%.16e\n", corr[tt * HM_N_CORR + k]);
                }
                tt = (t0 + T_glob / 2) % T_glob;
                fprintf(ofs, "%zu  %zu  %zu  %zu  %d  %.16e  %.16e\n", hi, hj, g1, g2, t,
                        corr[tt * HM_N_CORR + k], 0.0);
              }
            }
          }
        }
        fclose(ofs);
      }

      tm_stopwatch_pop(&g_timers, 0, 1, "");
    }  // for(max_time_slices)
  }  // for(max_samples)

  // freeing up memory
  for (size_t db = 0; db < 2; db++) {            // doublet: light or heavy
    for (size_t F = 0; F < 2; F++) {             // flavor dilution index
      for (size_t i_eo = 0; i_eo < 2; i_eo++) {  // even-odd index
        for (size_t i_f = 0; i_f < 2; i_f++) {   // flavor index of the doublet
          free(src_eo[db][F][i_eo][i_f]);
          free(prop_eo[db][F][i_eo][i_f]);
        }
      }
    }
  }

  free(phi1);
  free(phi2);
  free(corr_thr);
  free(corr);
#ifdef TM_USE_MPI
  free(corr_loc);
  free(corr_red);
#endif

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
