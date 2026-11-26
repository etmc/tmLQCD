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
#include "fatal_error.h"
#include "gamma.h"
#include "geometry_eo.h"
#include "gettime.h"
#include "global.h"
#include "invert_eo.h"
#include "io/params.h"
#include "io/spinor.h"
#include "io/utils.h"
#include "linalg/convert_eo_to_lexic.h"
#include "little_D.h"
#include "measurements.h"
#include "operator.h"
#include "operator/tm_operators_nd.h"
#include "ranlxs.h"
#include "solver/solver.h"
#include "source_generation.h"
#include "start.h"
#include "su3spinor.h"

#define TM_OMEAS_FILENAME_LENGTH 100

// whether the correlator measurement id uses operator i: 0 is the light quark, 1 the heavy
// doublet (HeavyMesons = yes), the operators in OSMesonsStrangeOperators the OS strange quarks
// (OSMesons = yes)
static bool correlators_uses_operator(const int id, const int i) {
  const measurement *const meas = &measurement_list[id];
  if (i == 0 || (i == 1 && meas->measure_heavy_mesons)) {
    return true;
  }
  for (int k = 0; meas->measure_os_mesons && k < meas->no_os_strange_operators; k++) {
    if (meas->os_strange_operators[k] == i) {
      return true;
    }
  }
  return false;
}

// checks the input of the Osterwalder-Seiler (OS) strange quarks of the correlator measurement id
// (OSMesons = yes), the type of the light quark operator 0 is checked before: fatal error if
// OSMesonsStrangeOperators is not given, if the light quark has 2KappaMu = 0, or if an entry is 0,
// not defined, not of type TMWILSON or CLOVER, has 2KappaMu = 0 (the default if 2KappaMu is
// missing) or 2KappaMu of the opposite sign of the light quark, or is listed more than once.
// At maximal twist, types 1 and 2 of the output (K+-, K0) are pseudoscalar correlators only if s
// is twisted like the light quark (u) and s' like d; otherwise (opposite sign, or an untwisted
// light quark) they are not (opposite sign: scalar correlators). With print_warnings, warnings if
// kappa or c_sw differ from the light quark.
static void os_mesons_check_input(const int id, const bool print_warnings) {
  const measurement *const meas = &measurement_list[id];
  char error[500];
  if (meas->no_os_strange_operators == 0) {
    snprintf(error, sizeof(error),
             "OSMesons = yes, but no OSMesonsStrangeOperators given in measurement %d", id);
    fatal_error(error, __func__);
  }
  const operator*const optr_l = &operator_list[0];
  if (optr_l->mu == 0.0) {
    snprintf(error, sizeof(error),
             "OSMesons = yes needs a twisted light quark: the light quark operator 0 has "
             "2KappaMu = 0 (measurement %d)",
             id);
    fatal_error(error, __func__);
  }
  // c_sw of the clover term (none for TMWILSON and WILSON)
  const double csw_l = (optr_l->type == CLOVER) ? optr_l->c_sw : 0.0;
  for (int k = 0; k < meas->no_os_strange_operators; k++) {
    const int i = meas->os_strange_operators[k];
    const char *problem = NULL;
    if (i == 0) {
      problem = "is the light quark operator";
    } else if (i >= no_operators) {
      problem = "is not defined in the input file";
    } else if (operator_list[i].type != TMWILSON && operator_list[i].type != CLOVER) {
      problem = "is not of type TMWILSON or CLOVER";
    } else if (operator_list[i].mu == 0.0) {
      problem = "has 2KappaMu = 0 (2KappaMu missing?)";
    } else if (operator_list[i].mu * optr_l->mu < 0.0) {
      problem =
          "has 2KappaMu of the opposite sign of the light quark operator 0 (the sign of s must be "
          "the one of the light quark, s' gets the opposite sign)";
    } else {
      for (int l = 0; l < k; l++) {
        if (meas->os_strange_operators[l] == i) {
          problem = "is listed more than once";
        }
      }
    }
    if (problem != NULL) {
      snprintf(error, sizeof(error),
               "OS strange quark operator %d (OSMesonsStrangeOperators of measurement %d) %s", i,
               id, problem);
      fatal_error(error, __func__);
    }

    const operator*const optr_s = &operator_list[i];
    const double csw_s = (optr_s->type == CLOVER) ? optr_s->c_sw : 0.0;
    if (print_warnings && g_proc_id == 0 &&
        (optr_s->kappa != optr_l->kappa || csw_s != csw_l)) {
      fprintf(stderr,
              "Warning! OS strange quark operator %d: kappa = %.16e, c_sw = %.16e differ from "
              "the light quark operator 0: kappa = %.16e, c_sw = %.16e (c_sw = 0: no clover "
              "term)!\n",
              i, optr_s->kappa, csw_s, optr_l->kappa, csw_l);
    }
  }
}

// checks the input of the correlator measurement id: fatal error if the light quark operator 0 is
// not defined or not of type TMWILSON, WILSON or CLOVER, if HeavyMesons = yes and the heavy
// doublet operator 1 is not defined or not of type DBTMWILSON or DBCLOVER, or for invalid OS
// strange quark input (OSMesons = yes, see os_mesons_check_input). With print_warnings, also the
// warnings of os_mesons_check_input and a warning for operators that are not used.
void correlators_check_input(const int id, const bool print_warnings) {
  const measurement *const meas = &measurement_list[id];
  char error[500];
  if (no_operators < 1 || (operator_list[0].type != TMWILSON && operator_list[0].type != WILSON &&
                           operator_list[0].type != CLOVER)) {
    snprintf(error, sizeof(error),
             "the correlator measurement %d needs the light quark operator 0 of type TMWILSON, "
             "WILSON or CLOVER%s",
             id, no_operators < 1 ? ", but no operator is defined in the input file" : "");
    fatal_error(error, __func__);
  }
  if (meas->measure_heavy_mesons && (no_operators < 2 || (operator_list[1].type != DBTMWILSON &&
                                                          operator_list[1].type != DBCLOVER))) {
    snprintf(error, sizeof(error),
             "HeavyMesons = yes needs the heavy doublet operator 1 of type DBTMWILSON or "
             "DBCLOVER%s (measurement %d)",
             no_operators < 2 ? ", but it is not defined in the input file" : "", id);
    fatal_error(error, __func__);
  }
  if (meas->measure_os_mesons) {
    os_mesons_check_input(id, print_warnings);
  }
  if (print_warnings && g_proc_id == 0) {
    bool first = true;
    for (int i = 1; i < no_operators; i++) {
      if (!correlators_uses_operator(id, i)) {
        fprintf(stderr, first ? "Warning! operator(s) %d" : ", %d", i);
        first = false;
      }
    }
    if (!first) {
      fprintf(stderr,
              " not used by the correlator online measurement (it uses operator 0 for the light "
              "quark, 1 for the heavy doublet with HeavyMesons = yes and OSMesonsStrangeOperators "
              "with OSMesons = yes)!\n");
    }
  }
}

static void os_mesons_measurement(const int id, const char *const suffix, const int t0,
                                  spinor *const xi[2], spinor *const psi_u[2]);

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
  // light quark operator 0 (type checked by correlators_check_input)
  optr = &operator_list[0];
  // we don't want to do inversion twice for this purpose here
  optr->DownProp = 0;

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
      // output filename (the suffix is also used for the OS-meson files)
      char suffix[32];
      if (max_samples == 1 && max_time_slices == 1) {
        snprintf(suffix, sizeof(suffix), "%06d", traj);
      } else if (max_samples == 1 && max_time_slices > 1) {
        snprintf(suffix, sizeof(suffix), "t%03d.%06d", ts, traj);
      } else {
        snprintf(suffix, sizeof(suffix), "s%03d.%06d", sample, traj);
      }
      snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "onlinemeas.%s", suffix);
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

      // OS mesons: same source and light propagator (the inversion above does not change the
      // source, the contraction above does not change the propagator)
      if (measurement_list[id].measure_os_mesons) {
        spinor *const xi[2] = {g_spinor_field[0], g_spinor_field[1]};
        spinor *const psi_u[2] = {g_spinor_field[2], g_spinor_field[3]};
        os_mesons_measurement(id, suffix, t0, xi, psi_u);
      }
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

  // dummy spinors for the heavy doublet inversion (allocated below)
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

  // selecting the operators i1 and i2 from the list of operators initialized before
  // (their types are checked by correlators_check_input)
  operator* optr1 = & operator_list[i1];  // light doublet
  operator* optr2 = & operator_list[i2];  // heavy c,s doublet
  // as in light_correlators_measurement: with AddDownPropagator = yes, op_invert would
  // additionally invert at -mu and return that propagator instead of the +mu one
  optr1->DownProp = 0;

  // allocating memory
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

// OS-meson correlators per time slice, k = 0..3 (output types 1..4): K+-, K0, connected eta_s,
// s-sbar', see os_mesons_measurement
#define OS_N_CORR 4

// contractions of the OS-meson measurement for one strange quark: adds the sums over the local
// spatial sites of every local time slice t to corr_thr[thread][t][k], i.e. per-thread partial
// sums (signs and normalisations are applied after the sums over the threads and ranks).
// Even-odd fields (index 0: even, 1: odd), see os_mesons_measurement:
// psi_u = D_u^{-1} xi, psi_s = D_s^{-1} xi, phi_s = D_s'^{-1} gamma_5 xi
static void os_mesons_contractions(double *const corr_thr, spinor *const psi_u[2],
                                   spinor *const psi_s[2], spinor *const phi_s[2]) {
#ifdef TM_USE_OMP
#pragma omp parallel
#endif
  {
#ifdef TM_USE_OMP
    double *const acc = corr_thr + (size_t)omp_get_thread_num() * T * OS_N_CORR;
#else
    double *const acc = corr_thr;
#endif
    for (int t = 0; t < T; t++) {
      const int j = g_ipt[t][0][0][0];  // source and sink separated by "t"
      double sum[OS_N_CORR] = {0.0};
      // static schedule: the same partition of the sites for a given number of threads
#ifdef TM_USE_OMP
#pragma omp for schedule(static) nowait
#endif
      for (int i = j; i < j + LX * LY * LZ; i++) {
        // parity and index of the lexicographic site i (as in convert_eo_to_lexic)
        const size_t eo = (g_lexic2eo[i] < (VOLUME + RAND) / 2) ? 0 : 1;
        const int ix = g_lexic2eosub[i];
        const spinor *const u = &psi_u[eo][ix];
        const spinor *const s = &psi_s[eo][ix];
        spinor g5_phi;
        _gamma5(g5_phi, phi_s[eo][ix]);
        // real parts of the scalar products; eta_s uses Re(phi^+ g5 psi_s) = Re(psi_s^+ g5 phi)
        sum[0] += _spinor_prod_re(*u, *s);      // K+-
        sum[1] += _spinor_prod_re(*u, g5_phi);  // K0
        sum[2] += _spinor_prod_re(*s, g5_phi);  // eta_s (connected)
        sum[3] += _spinor_prod_re(*s, *s);      // s-sbar'
      }
      for (size_t k = 0; k < OS_N_CORR; k++) {
        acc[t * OS_N_CORR + k] += sum[k];
      }
    }
  }  // omp parallel
}

/******************************************************
 *
 * Osterwalder-Seiler (OS) meson correlators, see the OS paper:
 * K. Ottnad et al., Nucl. Phys. B 896 (2015) 470, eqs. (4), (5), (7)
 *
 * One OS strange valence quark per operator in OSMesonsStrangeOperators (twisted mass
 * D_q = D_W + m_0 + i mu_q gamma_5): s with +mu_s and s' with -mu_s, where 2*kappa*mu_s is the
 * 2KappaMu of the operator (same sign as the light quark, see os_mesons_check_input). The source
 * xi (time slice t0, undiluted) and the light propagator psi_u = D_u^{-1} xi (light operator 0,
 * +mu_l) of the light measurement are reused (one-end trick, the d quark (-mu_l) enters via
 * gamma_5 hermiticity: D_d^{-1} = g5 (D_u^{-1})^+ g5).
 * Per strange operator two inversions:
 *   psi_s = D_s^{-1} xi,  phi_s = D_s'^{-1} gamma_5 xi
 * (normalised with 2*kappa as in op_invert). At maximal twist the pseudoscalar densities of the
 * physical basis are pseudoscalar densities in the twisted basis for (s, d) and (s, s'), and
 * scalar densities for (s', d) and (s, s), which gives the correlators (V_s: spatial volume, Re:
 * the imaginary parts are lattice artefacts and are dropped, see the OS paper)
 *   type 1 K+-   (sbar  i g5 d):  + Re sum_x psi_u^+ psi_s          / (V_s 2 kappa_l kappa_s)
 *   type 2 K0    (s'bar i g5 d):  - Re sum_x psi_u^+ g5 phi_s       / (V_s 2 kappa_l kappa_s)
 *   type 3 eta_s (sbar  i g5 s, connected only):
 *                                 - Re sum_x psi_s^+ g5 phi_s       / (V_s 2 kappa_s kappa_s)
 *   type 4 s-sbar' (sbar i g5 s'): + sum_x psi_s^+ psi_s            / (V_s 2 kappa_s kappa_s)
 * Type 1 has the normalisation of the light Cpp (onlinemeas type 1): with s = a copy of the
 * light operator, types 1 and 4 equal Cpp and types 2 and 3 are the connected neutral pion.
 * Output: os_mesons.opXX.<suffix> (XX: operator index) in the format of onlinemeas.<suffix>.
 *
 ******************************************************/
static void os_mesons_measurement(const int id, const char *const suffix, const int t0,
                                  spinor *const xi[2], spinor *const psi_u[2]) {
  tm_stopwatch_push(&g_timers, __func__, "");
  const measurement *const meas = &measurement_list[id];
  const double kappa_l = operator_list[0].kappa;
  char filename[TM_OMEAS_FILENAME_LENGTH];

  // source gamma_5 xi of s' and the strange propagators: passed to the inverters, whose
  // Hopping_Matrix exchanges the boundary into the field itself without halfspinor
  // -> they need the halo, i.e. VOLUMEPLUSRAND/2 (calloc: zero-initialised)
  spinor *g5_xi[2], *psi_s[2], *phi_s[2];
  for (int i_eo = 0; i_eo < 2; i_eo++) {
    g5_xi[i_eo] = calloc(VOLUMEPLUSRAND / 2, sizeof(spinor));
    psi_s[i_eo] = calloc(VOLUMEPLUSRAND / 2, sizeof(spinor));
    phi_s[i_eo] = calloc(VOLUMEPLUSRAND / 2, sizeof(spinor));
    gamma5(g5_xi[i_eo], xi[i_eo], VOLUME / 2);
  }
#ifdef TM_USE_OMP
  const int nthreads = omp_get_max_threads();
#else
  const int nthreads = 1;
#endif
  double *corr_thr = calloc((size_t)nthreads * T * OS_N_CORR, sizeof(double));
  double *corr = calloc((size_t)g_nproc_t * T * OS_N_CORR, sizeof(double));
#ifdef TM_USE_MPI
  double *corr_loc = calloc(T * OS_N_CORR, sizeof(double));  // sum over the threads
  double *corr_red = calloc(T * OS_N_CORR, sizeof(double));  // ... and over the ranks of t
#else
  double *corr_loc = corr;  // without MPI the sum over the threads is the result
#endif

  // source (time slice t0) and light propagator: one record each
  if (meas->write_os_mesons_sources) {
    tm_stopwatch_push(&g_timers, "writing_sources", "");
    snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "os_mesons_source.%s", suffix);
    write_eo_spinor_field(filename, xi[0], xi[1], 0);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
  }
  if (meas->write_os_mesons_props) {
    tm_stopwatch_push(&g_timers, "writing_propagators", "");
    snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "os_mesons_prop_light.%s", suffix);
    write_eo_spinor_field(filename, psi_u[0], psi_u[1], 0);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
  }

  for (int k = 0; k < meas->no_os_strange_operators; k++) {
    const int i_op = meas->os_strange_operators[k];
    operator* const optr = & operator_list[i_op];
    const double mu_s = optr->mu;
    if (g_debug_level > 1 && g_proc_id == 0) {
      printf("# OS strange quark operator %d: kappa = %.12f, mu = %.12f (s), %.12f (s')\n", i_op,
             optr->kappa, mu_s / 2. / optr->kappa, -mu_s / 2. / optr->kappa);
    }
    // the settings of the operator are restored after the inversions
    const int down_prop = optr->DownProp;
    spinor *const sr0 = optr->sr0, *const sr1 = optr->sr1;
    spinor *const prop0 = optr->prop0, *const prop1 = optr->prop1;
    // with AddDownPropagator = yes, op_invert would additionally invert at -mu and return that
    // propagator instead of the +mu one
    optr->DownProp = 0;

    tm_stopwatch_push(&g_timers, "inversions", "");
    // the propagator fields are the initial guess of the inversions (e.g. CG): start from zero
    for (int i_eo = 0; i_eo < 2; i_eo++) {
      zero_spinor_field(psi_s[i_eo], VOLUME / 2);
      zero_spinor_field(phi_s[i_eo], VOLUME / 2);
    }
    // s (+mu_s) on xi
    optr->sr0 = xi[0];
    optr->sr1 = xi[1];
    optr->prop0 = psi_s[0];
    optr->prop1 = psi_s[1];
    optr->inverter(i_op, 0, 0);
    // s' (-mu_s) on gamma_5 xi (as for AddDownPropagator = yes in op_invert, the deflation
    // subspace has to be updated when mu changes)
    optr->mu = -mu_s;
    dfl_subspace_updated = 1;
    optr->sr0 = g5_xi[0];
    optr->sr1 = g5_xi[1];
    optr->prop0 = phi_s[0];
    optr->prop1 = phi_s[1];
    optr->inverter(i_op, 0, 0);
    optr->mu = mu_s;
    dfl_subspace_updated = 1;
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    optr->DownProp = down_prop;
    optr->sr0 = sr0;
    optr->sr1 = sr1;
    optr->prop0 = prop0;
    optr->prop1 = prop1;

    // strange propagators: 2 records, psi_s (s on xi) and phi_s (s' on gamma_5 xi)
    if (meas->write_os_mesons_props) {
      tm_stopwatch_push(&g_timers, "writing_propagators", "");
      snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "os_mesons_prop_strange.op%02d.%s", i_op,
               suffix);
      write_eo_spinor_field(filename, psi_s[0], psi_s[1], 0);
      write_eo_spinor_field(filename, phi_s[0], phi_s[1], 1);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    }

    tm_stopwatch_push(&g_timers, "contractions", "");
    for (size_t l = 0; l < (size_t)nthreads * T * OS_N_CORR; l++) {
      corr_thr[l] = 0.0;
    }
    os_mesons_contractions(corr_thr, psi_u, psi_s, phi_s);
    // sum over the threads (in thread order)
    for (size_t l = 0; l < (size_t)T * OS_N_CORR; l++) {
      double sum = 0.0;
      for (int thr = 0; thr < nthreads; thr++) {
        sum += corr_thr[(size_t)thr * T * OS_N_CORR + l];
      }
      corr_loc[l] = sum;
    }
#ifdef TM_USE_MPI
    // sum over the ranks of each local time slice and gather of all time slices on the rank
    // that writes the output
    MPI_Reduce(corr_loc, corr_red, T * OS_N_CORR, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    if (g_mpi_time_rank == 0) {
      MPI_Gather(corr_red, T * OS_N_CORR, MPI_DOUBLE, corr, T * OS_N_CORR, MPI_DOUBLE, 0,
                 g_mpi_SV_slices);
    }
#endif
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    /* normalise and write everything into a file */
    if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
      const int T_glob = g_nproc_t * T;
      const double vol_fact = (g_nproc_x * LX) * (g_nproc_y * LY) * (g_nproc_z * LZ);
      const double kappa_s = optr->kappa;
      // sign and the kappa of the two quarks of each correlator (see above)
      const double sign[OS_N_CORR] = {+1.0, -1.0, -1.0, +1.0};
      const double kappa_q[OS_N_CORR] = {kappa_l, kappa_l, kappa_s, kappa_s};
      for (int t = 0; t < T_glob; t++) {
        for (size_t l = 0; l < OS_N_CORR; l++) {
          double *const c = &corr[(size_t)t * OS_N_CORR + l];
          *c = sign[l] * *c / vol_fact / 2.0 / kappa_q[l] / kappa_s;
        }
      }

      snprintf(filename, TM_OMEAS_FILENAME_LENGTH, "os_mesons.op%02d.%s", i_op, suffix);
      FILE *ofs = fopen(filename, "w");
      for (size_t l = 0; l < OS_N_CORR; l++) {
        const int type = (int)l + 1;
        int t, tt;
        fprintf(ofs, "%d  1  0  %.16e  %.16e\n", type, corr[t0 * OS_N_CORR + l], 0.0);
        for (t = 1; t < T_glob / 2; t++) {
          tt = (t0 + t) % T_glob;
          fprintf(ofs, "%d  1  %d  %.16e  ", type, t, corr[tt * OS_N_CORR + l]);
          tt = (t0 + T_glob - t) % T_glob;
          fprintf(ofs, "%.16e\n", corr[tt * OS_N_CORR + l]);
        }
        tt = (t0 + T_glob / 2) % T_glob;
        fprintf(ofs, "%d  1  %d  %.16e  %.16e\n", type, t, corr[tt * OS_N_CORR + l], 0.0);
      }
      fclose(ofs);
    }
  }  // strange operators

  for (int i_eo = 0; i_eo < 2; i_eo++) {
    free(g5_xi[i_eo]);
    free(psi_s[i_eo]);
    free(phi_s[i_eo]);
  }
  free(corr_thr);
  free(corr);
#ifdef TM_USE_MPI
  free(corr_loc);
  free(corr_red);
#endif
  tm_stopwatch_pop(&g_timers, 0, 1, "");
}

void correlators_measurement(const int traj, const int id, const int ieo) {
  // input checked at the start (init_measurements), again here without warnings since it may have
  // been changed in the meantime (hmc.reread)
  correlators_check_input(id, false);
  if (measurement_list[id].measure_heavy_mesons == 1) {
    const int i1 = 0, i2 = 1;  // light quark and heavy doublet operators
    heavy_correlators_measurement(traj, id, ieo, i1, i2);
  }

  light_correlators_measurement(traj, id, ieo);
}
