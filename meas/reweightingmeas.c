/***********************************************************************
 *
 * Measurements of the reweighting factors by Georg Bergner 2016
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
#include "measurements.h"
#include "correlators.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "linalg/mul_diff_mul_r.h"
#include "../operator/Hopping_Matrix.h"
#include "../operator/tm_operators.h"
#include "../operator/clovertm_operators.h"
#include "../operator/tm_operators_32.h"
#include "../operator/clovertm_operators_32.h"
#include "../operator/clover_leaf.h"
#include "../read_input.h"
#include "reweightingmeas.h"
#include "../operator.h"
#include "../DDalphaAMG_interface.h"
#include "../boundary.h"
#include "../global.h"
#include "solver/jdher.h"

#include <errno.h>
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif

/*#define CHECK_OPERATOR*/

typedef struct
{
  spinor* buffer;
  spinor** ar;
  unsigned int length;
} spinor_array;

static void
set_zero (spinor * const R, const int N)
{
#ifdef TM_USE_OMP
#pragma omp parallel
    {
#endif

  int ix;
  spinor *r;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for (ix = 0; ix < N; ++ix)
    {
      r = (spinor *) R + ix;

      r->s0.c0 = 0.0;
      r->s0.c1 = 0.0;
      r->s0.c2 = 0.0;

      r->s1.c0 = 0.0;
      r->s1.c1 = 0.0;
      r->s1.c2 = 0.0;

      r->s2.c0 = 0.0;
      r->s2.c1 = 0.0;
      r->s2.c2 = 0.0;

      r->s3.c0 = 0.0;
      r->s3.c1 = 0.0;
      r->s3.c2 = 0.0;
    }
#ifdef TM_USE_OMP
} /* OpenMP closing brace */
#endif
}

static int
alloc_spinor_field (spinor_array* ar, int evenodd)
{
  int i = 0;
  unsigned long int volume = VOLUMEPLUSRAND / 2;
  if (evenodd == 0)
    {
      volume = VOLUMEPLUSRAND;
    }
  if (ar->buffer != NULL || ar->ar != NULL || ar->length == 0)
    {
      return (3);
    }
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  if((void*)(ar->buffer = (spinor*)shmalloc((ar->length*volume+1)*sizeof(spinor))) == NULL)
    {
      printf ("malloc errno : %d\n",errno);
      errno = 0;
      return(1);
    }
#else
  if ((void*) (ar->buffer = (spinor*) calloc (ar->length * volume + 1,
                                              sizeof(spinor))) == NULL)
    {
      printf ("malloc errno : %d\n", errno);
      errno = 0;
      return (1);
    }
#endif
  if ((void*) (ar->ar = (spinor**) malloc (ar->length * sizeof(spinor*)))
      == NULL)
    {
      printf ("malloc errno : %d\n", errno);
      errno = 0;
      return (2);
    }
#if ( defined SSE || defined SSE2 || defined SSE3)
  ar->ar = (spinor*)(((unsigned long int)(ar->buffer)+ALIGN_BASE)&~ALIGN_BASE);
#else
  ar->ar[0] = ar->buffer;
#endif

  for (i = 1; i < ar->length; i++)
    {
      ar->ar[i] = g_spinor_field[i - 1] + volume;
    }

  return (0);
}

static void
free_spinor_field (spinor_array* ar)
{
  if (ar->buffer != NULL)
    {
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
      shfree(ar->buffer);
#else
      free (ar->buffer);
#endif
    }
  if (ar->ar != NULL)
    {
      free (ar->ar);
    }
  ar->ar = NULL;
  ar->buffer = NULL;
  ar->length = 0;
}

static int
is_schur_complement (const matrix_mult f)
{
  if (f == Msw_psi ||     //          Schur complement with mu=0 on odd sites
      f == Qsw_psi || // Gamma5 - Schur complement with mu=0 on odd sites
      f == Mtm_plus_psi ||  //          Schur complement with plus mu
      f == Msw_plus_psi ||  //          Schur complement with plus mu
      f == Qtm_plus_psi ||  // Gamma5 - Schur complement with plus mu
      f == Qsw_plus_psi ||  // Gamma5 - Schur complement with plus mu
      f == Mtm_minus_psi || //          Schur complement with minus mu
      f == Msw_minus_psi || //          Schur complement with minus mu
      f == Qtm_minus_psi || // Gamma5 - Schur complement with minus mu
      f == Qsw_minus_psi || // Gamma5 - Schur complement with minus mu
      f == Qtm_pm_psi ||    //          Schur complement squared
      f == Qsw_pm_psi)
    {   //          Schur complement squared
      return 1;
    }
  return 0;
}

static double
get_sw_reweighting (const double mu1, const double mu2, const double kappa1,
                    const double kappa2, const double csw1, const double csw2)
{
  double ret;
  sw_term ((const su3**) g_gauge_field, kappa1, csw1);
  ret = -sw_trace (0, mu1);
  if (kappa1 != kappa2 || csw1 != csw2)
    {
      sw_term ((const su3**) g_gauge_field, kappa2, csw2);
    }
  ret += sw_trace (0, mu2);
  return (ret);
}

static int
is_sym_pos_definite (const matrix_mult f)
{
  if (f == Qtm_pm_psi ||    //          Schur complement squared
      f == Qsw_pm_psi || f == Q_pm_psi)
    { //          Full operator    squared
      return 1;
    }
  return 0;
}

static void
update_global_parameters (const int op_id)
{
  operator * optr = &operator_list[op_id];
  g_kappa = optr->kappa;
  boundary (g_kappa);
  g_mu = optr->mu;
  g_c_sw = optr->c_sw;
  if (optr->type == CLOVER)
    {
      if (g_cart_id == 0 && g_debug_level > 1)
        {
          printf ("#\n# csw = %e, computing clover leafs\n", g_c_sw);
        }
      init_sw_fields (VOLUME);

      sw_term ((const su3**) g_gauge_field, optr->kappa, optr->c_sw);
      /* this must be EE here!   */
      /* to match clover_inv in Qsw_psi */
      sw_invert (EE, optr->mu);
      /* now copy double sw and sw_inv fields to 32bit versions */
      copy_32_sw_fields ();
    }/*clover leave update*/
}


static int
invert_operator_Q (spinor * const P, spinor * const Q, const int op_id,
                   const int pm, int updateparameters)
{

  operator * optr = &operator_list[op_id];
  int iteration_count = 0;
  int use_solver = optr->solver;
  int rel_prec = optr->rel_prec;
  double eps_sq = optr->eps_sq;
  double check_prec;
  int max_iter = optr->maxiter;
  matrix_mult f = optr->applyQm;
  int is_squared = 0;
  spinor * source = Q;
  if (pm == 1)
    {
      f = optr->applyQp;
    }
  solver_params_t solver_params = optr->solver_params;
  int n = VOLUME;
  if (is_schur_complement (f))
    {
      n = VOLUME / 2;
    }

  optr->iterations = 0;
  optr->reached_prec = -1.;
  update_global_parameters (op_id);

  if (use_solver == MIXEDCG || use_solver == RGMIXEDCG || use_solver == CG)
    {
      if (!is_sym_pos_definite (f))
        {
          f = optr->applyQsq;
          is_squared = 1;
          if (pm == 0)
            {
              optr->applyQp (g_spinor_field[DUM_DERI], Q);
              source = g_spinor_field[DUM_DERI];
            }
        }
    }

  /*check initial precision since some inverters fail if already converged*/
  f (g_spinor_field[DUM_DERI + 2], P);
  diff (g_spinor_field[DUM_DERI + 2], g_spinor_field[DUM_DERI + 2], source, n);
  check_prec = square_norm (g_spinor_field[DUM_DERI + 2], VOLUME, 1);
  if (g_proc_id == 0)
    {
      printf ("Inversion initial precision %e\n", check_prec);
    }
  if (check_prec < 1e-24)
    {
      return (0);
    }

  if (use_solver == MIXEDCG || use_solver == RGMIXEDCG)
    {
      // the default mixed solver is rg_mixed_cg_her
      int
      (*msolver_fp) (spinor * const, spinor * const, solver_params_t, const int,
                     double, const int, const int, matrix_mult,
                     matrix_mult32) = rg_mixed_cg_her;

      // but it might be necessary at some point to use the old version
      if (use_solver == MIXEDCG)
        {
          msolver_fp = mixed_cg_her;
        }

      if (usegpu_flag)
        {
#ifdef HAVE_GPU
#ifdef TEMPORALGAUGE
          to_temporalgauge(g_gauge_field, source , P);
#endif
          iteration_count = linsolve_eo_gpu(P, source, max_iter, eps_sq, rel_prec, n, f);
#ifdef TEMPORALGAUGE
          from_temporalgauge(source, P);
#endif
#endif
          return (iteration_count);
        }
      else
        {
          if (f == Qtm_pm_psi)
            {
              iteration_count = msolver_fp (P, source, solver_params, max_iter,
                                            eps_sq, rel_prec, n, f,
                                            &Qtm_pm_psi_32);
              return (iteration_count);
            }
          else if (f == Q_pm_psi)
            {
              iteration_count = msolver_fp (P, source, solver_params, max_iter,
                                            eps_sq, rel_prec, n, f,
                                            &Q_pm_psi_32);
              return (iteration_count);
            }
          else if (f == Qsw_pm_psi)
            {
              copy_32_sw_fields ();
              iteration_count = msolver_fp (P, source, solver_params, max_iter,
                                            eps_sq, rel_prec, n, f,
                                            &Qsw_pm_psi_32);
              return (iteration_count);
            }
          else
            {
              if (g_proc_id == 0)
                printf (
                    "Warning: 32 bit matrix not available. Falling back to CG in 64 bit\n");
              use_solver = CG;
            }
        }
    }
  if (use_solver == CG)
    {
      iteration_count = cg_her (P, source, max_iter, eps_sq, rel_prec, n, f);
    }
  else if (use_solver == BICGSTAB)
    {
      iteration_count = bicgstab_complex (P, source, max_iter, eps_sq, rel_prec,
                                          n, f);
    }
#ifdef DDalphaAMG
  else if (use_solver == MG)
    {
      if (updateparameters == 1)
        {
          MG_update_kappa (g_kappa);
        }
      iteration_count = MG_solver (P, source, eps_sq, max_iter, rel_prec, n,
                                   g_gauge_field, f);
    }
#endif
  else
    {
      if (g_proc_id == 0)
        printf (
            "Error: solver not allowed for degenerate solve. Aborting...\n");
      return -1;
    }
  if (is_squared)
    {
      f = optr->applyQm;
      if (pm == 1)
        {
          optr->applyQm (g_spinor_field[DUM_DERI], P);
          assign (P, g_spinor_field[DUM_DERI], n);
          f = optr->applyQp;
        }
    }
  /*check precision*/
  f (g_spinor_field[DUM_DERI], P);
  diff (g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], Q, n);
  check_prec = square_norm (g_spinor_field[DUM_DERI], VOLUME, 1);
  optr->reached_prec = check_prec;
  optr->iterations = iteration_count;
  if (g_proc_id == 0)
    {
      printf ("Inversion final precision %e\n", check_prec);
    }
  return (iteration_count);
}

//#define DEBUG_PARAMETER_CHANGE

static void
interpolate (double * const rcurrent, const double rinitial,
             const double rfinal, const int numsteps, const int thisstep)
{
#ifdef DEBUG_PARAMETER_CHANGE
  if(thisstep==1){
      (*rcurrent)=100;
  }else{
      (*rcurrent)=10;
  }
  return;
#endif
  (*rcurrent) = rinitial
      + (thisstep + 1) * (rfinal - rinitial) / ((double) numsteps);
}

static void
chebyshev_coeff (const unsigned int np, double * const coeff, const double a,
                 const double b)
{
  double * fxk;
  double * acxk;
  fxk = malloc (np * sizeof(double));
  acxk = malloc (np * sizeof(double));
  int n, k;
  double fakt;
  for (n = 0; n < np; n++)
    {
      acxk[n] = (3.14159265358979323846) * ((double) n + 0.5) / ((double) np);
      fxk[n] = log (a * cos (acxk[n]) + b);
    }
  for (k = 0; k < np; k++)
    {
      coeff[k] = 0;
      fakt = (k == 0) ? 1.0 / (double) np : 2.0 / (double) np;
      for (n = 0; n < np; n++)
        {
          coeff[k] += fakt * cos (k * acxk[n]) * fxk[n];
        }
    }
  free (fxk);
  free (acxk);
}

static double
poly_cheb (const unsigned int np, const double * const coeff, const double x)
{
  double y, t1, t0, t2;
  int j;
  y = coeff[0];
  if (np < 1)
    {
      return y;
    }
  t0 = 1.0;
  t1 = x;
  y += coeff[1] * t1;
  for (j = 1; j + 1 < np; j++)
    {
      t2 = 2.0 * x * t1 - t0;
      t0 = t1;
      t1 = t2;
      y += coeff[j + 1] * t1;
    }
  return y;
}

/*Trivial test tests the function with the scaled idenetity matrix.*/
/*#define TRIVIAL_TEST*/
/*Convergence chesk: comparison of order n and n+1.*/
#define CHECK_CONVERGENCE
static void
log_determinant_estimate (const int operatorid, int chebmax, int estimators,
                          const double minev, const double maxev,
                          const double kappa1, const double kappa2,
                          const double kappa2Mu1, const double kappa2Mu2,
                          const double shift, const int traj,
                          const split_list * const split,
                          const vector_list * const coefflist)
{
  double * coeff;
  const double t1 = maxev - minev;
  const double t2 = maxev + minev;
  const double a = 2.0 / t1;
  const double b = -t2 / t1;
  const double am = t1 / 2.0;
  const double bm = t2 / 2.0;
  double rel_shift;
  int k, l, n;
  int orderstart, orderend;
  int splitlength, sl;
  spinor_array spinorarray;
  static int written = 0;
  double x, y, y1, prodre;
  FILE * ofs;
  spinor * vs;
  spinor * u;
#ifdef CHECK_CONVERGENCE
  spinor * unm1;
  double prodrenm1;
#endif
  spinor * v0;
  spinor * v1;
  spinor * v2;
  spinor * vt0;
  spinor * vt1;
  spinor * vt2;
  spinor * tmp;
  operator * optr;
  char* filename;
  char buf[100];

  spinorarray.ar = NULL;
  spinorarray.buffer = NULL;
  spinorarray.length = 0;

  n = VOLUME;

  filename = buf;
  sprintf (filename, "%s%.6d", "reweightingmeas_cheb.", traj);

  optr = &operator_list[operatorid];

  spinorarray.length = 9;
  if (is_schur_complement (optr->applyQsq))
    {
      n = VOLUME / 2;
      alloc_spinor_field (&spinorarray, 1);
    }
  else
    {
      alloc_spinor_field (&spinorarray, 0);
    }

  vs = spinorarray.ar[0];
  u = spinorarray.ar[1];
  v0 = spinorarray.ar[2];
  v1 = spinorarray.ar[3];
  v2 = spinorarray.ar[4];
  vt0 = spinorarray.ar[5];
  vt1 = spinorarray.ar[6];
  vt2 = spinorarray.ar[7];
#ifdef CHECK_CONVERGENCE
  unm1 = spinorarray.ar[8];
#endif

  if (coefflist->el != NULL && coefflist->s != 0)
    {
      chebmax = coefflist->s;
      coeff = malloc (chebmax * sizeof(double));
      for (k = 0; k < chebmax; k++)
        {
          coeff[k] = coefflist->el[k];
        }
    }
  else
    {
      coeff = malloc (chebmax * sizeof(double));
      chebyshev_coeff (chebmax, coeff, am, bm);
    }
  if (g_proc_id == 0 && g_debug_level > 3 && written == 0)
    {
      written = 1;
      ofs = fopen ("polynomialapproxtest.txt", "w");
      fprintf (
          ofs,
          "# Test of the Chebyshev approximation of the log(x) function: index x y |y-log(x)| chebyshevorder minx, maxx\n");
      for (k = 0; k < 200; k++)
        {
          x = minev + (maxev - minev) * (double) k / (double) (200 - 1);
          y = poly_cheb (chebmax, coeff, a * x + b);
          y1 = NAN;
          if (chebmax > 0)
            {
              y1 = poly_cheb (chebmax - 1, coeff, a * x + b);
            }
          fprintf (ofs, "%d %g   %g %g   %g %g   %d %g %g\n", k, x, y, y1,
                   fabs (y - log (x)), fabs (y1 - log (x)), chebmax, minev,
                   maxev);
        }
      fclose (ofs);
      ofs = fopen ("coeff_out.txt", "w");
      for (k = 0; k < chebmax; k++)
        {
          fprintf (ofs, "%g\n", coeff[k]);
        }
      fclose (ofs);
    }

  /* include T_min(x) to T_(max-1) (x) */
  orderstart = 0;
  orderend = chebmax;
  splitlength = 1;
  rel_shift = shift;
  if (split->est != NULL && split->ord != NULL && split->s != 0)
    {
      splitlength = split->s;
      orderend = split->ord[0];
      estimators = split->est[0];
    }
  for (sl = 0; sl < splitlength; sl++)
    {
      if (sl > 0)
        {
          orderstart = orderend;
          orderend = split->ord[sl];
          estimators = split->est[sl];
          rel_shift = 0.0;
        }

      if (orderend > 0)
        {
          for (k = 0; k < estimators; k++)
            {
              if (orderstart > 1 || orderend < 2)
                {
                  set_zero (u, n);
#ifdef CHECK_CONVERGENCE
                  set_zero (unm1, n);
#endif
                }
              /*
               * Generate estimator
               */
              if (n == VOLUME / 2)
                {
                  random_spinor_field_eo (vs, reproduce_randomnumber_flag,
                                          RN_Z2);
                }
              else
                {
                  random_spinor_field_lexic (vs, reproduce_randomnumber_flag,
                                             RN_Z2);
                }

              assign (v0, vs, n);
              assign (vt0, vs, n);

#ifdef TRIVIAL_TEST
              prodre = scalar_prod_r(vs1, vs1, n, 1);
              printf("Test noise <v|v>/Volume= %g, Volume=%d\n",prodre/(double)n,n);
              mul_r(v1,kappa1,vs1,n);
#else
              optr->kappa = kappa1;
              optr->mu = kappa2Mu1;
              update_global_parameters (operatorid);
              optr->applyQsq (v1, vs);
#endif
              /* Makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
              assign_mul_add_mul_r (v1, vs, a, b, n);

#ifdef TRIVIAL_TEST
              mul_r(vt1,kappa2,vs1,n);
#else
              optr->kappa = kappa2;
              optr->mu = kappa2Mu2;
              update_global_parameters (operatorid);
              optr->applyQsq (vt1, vs);
#endif
              /* Makes (*R)=c1*(*R)+c2*(*S) , c1 and c2 are real constants */
              assign_mul_add_mul_r (vt1, vs, a, b, n);

              if (orderstart < 2)
                {
                  /* Makes (*R)=c1*(*S)-c2*(*U) , c1 and c2 are complex constants */
                  mul_diff_mul_r (u, vt1, v1, coeff[1], coeff[1], n);

#ifdef CHECK_CONVERGENCE
                  if (orderend > 1)
                    {
                      mul_diff_mul_r (unm1, vt1, v1, coeff[1], coeff[1], n);
                    }
#endif
                }

              /*This part is only more efficient if not the clover term h
               * as to be recompted everytime the parameters are changed.*/
#if 0
              for (l = 1; l + 1 < orderend; l++)
                {
#ifndef TRIVIAL_TEST
                  optr->kappa = kappa1;
                  optr->mu = kappa2Mu1;
                  update_global_parameters(operatorid);

                  optr->applyQsq(v2, v1);
#else
                  mul_r(v2,kappa1,v1,n);
#endif
                  if(l+1>=orderstart)
                    {
                      /* Makes (*R) = c1*(*R) + c2*(*S) + c3*(*U) */
                      assign_mul_add_mul_add_mul_r(v2, v1, v0, 2.0 * a, 2.0 * b, -1.0,
                          n);
                    }
                  tmp = v0;
                  v0 = v1;
                  v1 = v2;
                  v2 = tmp;

#ifndef TRIVIAL_TEST
                  optr->kappa = kappa2;
                  optr->mu = kappa2Mu2;
                  update_global_parameters(operatorid);

                  optr->applyQsq(vt2, vt1);
#else
                  mul_r(vt2,kappa2,vt1,n);
#endif
                  /* Makes (*R) = c1*(*R) + c2*(*S) + c3*(*U) */
                  assign_mul_add_mul_add_mul_r(vt2, vt1, vt0, 2.0 * a, 2.0 * b,
                      -1.0, n);
                  tmp = vt0;
                  vt0 = vt1;
                  vt1 = vt2;
                  vt2 = tmp;

                  if(l+1>=orderstart)
                    {
                      /* (*R) = (*R) + c1*(*S) + c2*(*U) */
                      assign_add_mul_add_mul_r(u, vt1, v1, coeff[l + 1],
                          -coeff[l + 1], n);
                    }

                }

#endif

#ifndef TRIVIAL_TEST
              optr->kappa = kappa1;
              optr->mu = kappa2Mu1;
              update_global_parameters (operatorid);
#endif

              for (l = 1; l + 1 < orderend; l++)
                {
#ifdef TRIVIAL_TEST
                  mul_r(v2,kappa1,v1,n);
#else
                  optr->applyQsq (v2, v1);
#endif

                  /* Makes (*R) = c1*(*R) + c2*(*S) + c3*(*U) */
                  assign_mul_add_mul_add_mul_r (v2, v1, v0, 2.0 * a, 2.0 * b,
                                                -1.0, n);
                  tmp = v0;
                  v0 = v1;
                  v1 = v2;
                  v2 = tmp;

                  if (l + 1 >= orderstart)
                    {
                      /*   (*P) = (*P) + c(*Q)        c is a real constant   */
                      assign_add_mul_r (u, v1, -coeff[l + 1], n);
#ifdef CHECK_CONVERGENCE
                      if (l + 2 < orderend)
                        {
                          assign_add_mul_r (unm1, v1, -coeff[l + 1], n);
                        }
#endif
                    }
                }
#ifndef TRIVIAL_TEST
              optr->kappa = kappa2;
              optr->mu = kappa2Mu2;
              update_global_parameters (operatorid);
#endif

              for (l = 1; l + 1 < orderend; l++)
                {
#ifdef TRIVIAL_TEST
                  mul_r(vt2,kappa2,vt1,n);
#else
                  optr->applyQsq (vt2, vt1);
#endif
                  /* Makes (*R) = c1*(*R) + c2*(*S) + c3*(*U) */
                  assign_mul_add_mul_add_mul_r (vt2, vt1, vt0, 2.0 * a, 2.0 * b,
                                                -1.0, n);
                  tmp = vt0;
                  vt0 = vt1;
                  vt1 = vt2;
                  vt2 = tmp;
                  if (l + 1 >= orderstart)
                    {
                      /*   (*P) = (*P) + c(*Q)        c is a real constant   */
                      assign_add_mul_r (u, vt1, coeff[l + 1], n);
#ifdef CHECK_CONVERGENCE
                      if (l + 2 < orderend)
                        {
                          assign_add_mul_r (unm1, vt1, coeff[l + 1], n);
                        }
#endif
                    }
                }

              prodre = scalar_prod_r (vs, u, n, 1);
#ifdef CHECK_CONVERGENCE
              prodrenm1 = scalar_prod_r (vs, unm1, n, 1);
#endif

              if (g_proc_id == 0)
                {
                  ofs = fopen (filename, "a");
                  fprintf (ofs, "%d %d %g %g %g %g %d %d %g %g ", traj, sl,
                           kappa1, kappa2, kappa2Mu1, kappa2Mu2, orderend,
                           orderstart, minev, maxev);
                  fprintf (ofs, "     %.17g %.17g ", prodre + rel_shift,
                           prodre);
#ifdef CHECK_CONVERGENCE
                  fprintf (ofs, "      %.17g \n", prodrenm1 + rel_shift);
#endif
                  fclose (ofs);
                }

            } /* estimator iteration*/
        } /* orderend>0*/
    } /* splitlength loop */
  free (coeff);
  free_spinor_field (&spinorarray);
}

void
reweighting_measurement (const int traj, const int id, const int ieo)
{
  reweighting_parameter* param;
  double atime, etime;
  operator * optr;
  FILE *ofs;
  FILE *ofs_full;
  char *filename;
  char *filename_full;
  char buf[100];
  char buf_full[100];
  double square1, square2, prodre, prodim, prodreg5, prodimg5, cswpart;
  double kappa0, kappa, k2mu0, k2mu, csw0, csw;
  double kappafinal, k2mufinal, cswfinal, rmufinal;
  double kappainitial, k2muinitial, cswinitial, rmuinitial;
  double mdiff, rmu, rmu0;
  double tmp1, tmp2, tmp3;
  int updateinverter;
  int kapparew;
  int murew;
#ifdef CHECK_OPERATOR
  double checkg5,check1,check2,check3,check4,check5;
  double kappa_old;
#endif
  /* now we bring it to normal format */
  /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
  spinor * lexicfield1;
  spinor * lexicfield2;
#ifdef CHECK_OPERATOR
  spinor * lexicfield3;
  spinor * lexicfield4;
#endif
  int operatorid;
  int numsamples, snum, internum;
  unsigned long int site;

  updateinverter=0;
  param = (reweighting_parameter*) measurement_list[id].parameter;
  lexicfield1 = g_spinor_field[4];
  lexicfield2 = g_spinor_field[6];
#ifdef CHECK_OPERATOR
  lexicfield3=g_spinor_field[8];
  lexicfield4=g_spinor_field[10];
#endif
  operatorid = param->reweighting_operator;
  numsamples = param->reweighting_number_sources;
  filename = buf;
  filename_full = buf_full;
  sprintf (filename, "%s%.6d", "reweightingmeas.", traj);
  sprintf (filename_full, "%s%.6d", "reweightingmeas_full_data.", traj);
  if (g_proc_id == 0)
    {
      fprintf (stdout, "Reweighting measurement %d with %d samples.\n", id,
               numsamples);
    }

  init_operators ();
  if (no_operators < operatorid)
    {
      if (g_proc_id == 0)
        {
          fprintf (
              stderr,
              "Warning! Number of operators smaller than the given number for the reweighting operator, unable to perform measurement!\n");
        }
      return;
    }
  atime = gettime ();

  optr = &operator_list[operatorid];
  optr->DownProp = 0;
  optr->sr0 = g_spinor_field[0];
  optr->sr1 = g_spinor_field[1];
  optr->prop0 = g_spinor_field[2];
  optr->prop1 = g_spinor_field[3];

  /* now checking parameters
   * The target values are obtained from the operator*/

  kappafinal = optr->kappa;
  k2mufinal = optr->mu;
  cswfinal = optr->c_sw;

  /* the initial values are obtained
   * from the parameter
   */
  k2muinitial = param->k2mu0;
  kappainitial = param->kappa0;
  cswinitial = 0;
  kapparew = 1;
  if (kappainitial == kappafinal)
    {
      kapparew = 0;
    }
  if (kappainitial == 0.0)
    {
      kappainitial = kappafinal;
      kapparew = 0;
    }
  if (cswinitial == 0.0)
    {
      cswinitial = cswfinal;
    }
  /* be careful:
   * in the mu reweighting it is the parameter Mu and not 2KappaMu that
   * counts
   */
  rmufinal = k2mufinal / 2.0 / kappafinal;
  rmuinitial = k2muinitial / 2.0 / kappainitial;
  murew = 1;
  if (k2muinitial == 0.0)
    {
      rmuinitial = rmufinal;
      k2muinitial = 2.0 * kappainitial * rmuinitial;
      murew = 0;
    }
  if (fabs (rmuinitial - rmufinal) < 1e-14)
    {
      murew = 0;
    }

  /* second option:
   * determine mu and mu0 explicitly in the
   */
  if (g_proc_id == 0 && (param->rmu != 0 || param->rmu0 != 0))
    {
      printf (
          "WARNING: measurement uses custom mu values mu=%e, mu0=%e instead of the parameters of the operator.\n",
          param->rmu, param->rmu0);
    }
  if (param->rmu0 != 0)
    {
      rmuinitial = param->rmu0;
      k2muinitial = 2.0 * kappainitial * rmuinitial;
    }
  if (param->rmu != 0)
    {
      rmufinal = param->rmu;
      k2mufinal = 2.0 * kappafinal * rmufinal;
    }
  if (fabs (rmuinitial - rmufinal) > 1e-14)
    {
      murew = 1;
    }

  if (murew && (g_proc_id == 0))
    {
      printf ("Mu reweighting chosen: ");
      printf ("mu=%e to mu=%e.\n", rmuinitial, rmufinal);
    }

  if (kapparew && (g_proc_id == 0))
    {
      printf ("Kappa reweighting chosen: ");
      printf ("kappa=%e to kappa=%e\n", kappainitial, kappafinal);
    }

  if (!murew && !kapparew)
    {
      if (g_proc_id == 0)
        {
          printf ("ERROR: no mu or kappa reweighting.\n");
        }
    }

  if (murew && kapparew)
    {
      if (g_proc_id == 0)
        {
          printf (
              "WARNING: Mu AND Kappa reweighting chosen. If unprecond version, the rawdata has to be combined by hand!\n");
        }
    }

  if (param->use_evenodd && !is_schur_complement (optr->applyQsq)
      && g_proc_id == 0)
    {
      printf (
          "WARNING: If you want to use the preconditioned version the operator should be even-odd preconditioned.\n");
    }

  cswpart = 0;

  if (param->evest)
    {
      if (g_proc_id == 0)
        {
          printf ("Calculating minimal/maximal eigenvalue estimates -- no longer supported!\n");
        }
      //estimate_eigenvalues (operatorid, traj);
    }

  if (param->use_cheb)
    {
      if (is_schur_complement (optr->applyQsq) && !param->use_evenodd
          && g_proc_id == 0)
        {
          printf (
              "WARNING: If you want to use Chebyshev approximation without preconditioning you have to switch it of in the operator as well.\n");
        }
      if (param->use_evenodd == 1)
        {
          cswpart = get_sw_reweighting (k2muinitial, k2mufinal, kappainitial,
                                        kappafinal, cswinitial, cswfinal);
        }
      log_determinant_estimate (operatorid, param->cheborder,
                                param->estimatorscheb, param->minev,
                                param->maxev, kappainitial, kappafinal,
                                k2muinitial, k2mufinal, cswpart, traj,
                                &param->splitlist, &param->coeff);
      if (param->only_cheb)
        {
          return;
        }
    }

  kappa0 = kappa = kappainitial;
  k2mu0 = k2mu = k2muinitial;
  csw0 = csw = cswinitial;
  rmu0 = rmu = rmuinitial;

  if (param->interpolationsteps < 1)
    param->interpolationsteps = 1;

  for (internum = 0; internum < param->interpolationsteps; internum++)
    {
      if (kapparew)
        {
          kappa0 = kappa;
          tmp1 = 1.0 / kappainitial;
          tmp2 = 1.0 / kappafinal;
          tmp3 = 1.0 / kappa;
          interpolate (&tmp3, tmp1, tmp2, param->interpolationsteps, internum);
          kappa = 1.0 / tmp3;
          updateinverter=1;
        }
      if (murew)
        {
          /* use quadratic interpolation in the case of mu*/
          rmu0 = rmu;
          tmp1 = rmuinitial * rmuinitial;
          tmp2 = rmufinal * rmufinal;
          tmp3 = rmu * rmu;
          interpolate (&tmp3, tmp1, tmp2, param->interpolationsteps, internum);
          rmu = sqrt (tmp3);
          updateinverter=1;
        }
      k2mu = 2.0 * kappa * rmu;
      k2mu0 = 2.0 * kappa0 * rmu0;
      optr->kappa = kappa;
      optr->mu = k2mu;
      optr->c_sw = csw;

      if (param->use_evenodd == 1)
        {
          cswpart = get_sw_reweighting (k2mu0, k2mu, kappa0, kappa, csw0, csw);
        }

      for (snum = 0; snum < numsamples; snum++)
        {
          if (param->use_evenodd == 0)
            {
              random_spinor_field_eo (optr->sr0, reproduce_randomnumber_flag,
                                      RN_GAUSS);
              random_spinor_field_eo (optr->sr1, reproduce_randomnumber_flag,
                                      RN_GAUSS);
              convert_eo_to_lexic (lexicfield1, optr->sr0, optr->sr1);
              square1 = square_norm (lexicfield1, VOLUME, 1);
              // we don't want to do inversion twice for this purpose here
              // op_id = 0, index_start = 0, write_prop = 0
              optr->inverter (operatorid, 0, 0);
              convert_eo_to_lexic (lexicfield2, optr->prop0, optr->prop1);
              square2 = square_norm (lexicfield2, VOLUME, 1);
              prodre = scalar_prod_r (lexicfield1, lexicfield2, VOLUME, 1);
              prodim = scalar_prod_i (lexicfield1, lexicfield2, VOLUME, 1);
              for (site = 0; site < VOLUME; site++)
                {
                  _gamma5(lexicfield2[site], lexicfield2[site]);
                }
              prodreg5 = scalar_prod_r (lexicfield1, lexicfield2, VOLUME, 1);
              prodimg5 = scalar_prod_i (lexicfield1, lexicfield2, VOLUME, 1);
#ifdef CHECK_OPERATOR
              for(site = 0; site < VOLUME; site++)
                {
                  _gamma5(lexicfield2[site], lexicfield1[site]);
                }
              checkg5=scalar_prod_r(lexicfield1,lexicfield2,VOLUME, 1); /* should be zero*/

              kappa_old=optr->kappa;
              optr->applyM(optr->prop0, optr->prop1, optr->sr0, optr->sr1); /* prop=D rand*/
              mul_r(optr->prop0, (0.5/optr->kappa), optr->prop0, VOLUME / 2);/*correct normalisation*/
              mul_r(optr->prop1, (0.5/optr->kappa), optr->prop1, VOLUME / 2);
              convert_eo_to_lexic(lexicfield2, optr->prop0, optr->prop1);
              check1 = -scalar_prod_r(lexicfield1,lexicfield2,VOLUME, 1);
              check2 = -scalar_prod_i(lexicfield1,lexicfield2,VOLUME, 1);
              check3 = scalar_prod_r(lexicfield2,lexicfield2,VOLUME ,1);

              optr->kappa=1.0/(1.0/kappa_old+2.0);
              optr->inverter(operatorid,0,0); /* to ensure that everything is updated*/
              optr->applyM(optr->prop0, optr->prop1, optr->sr0, optr->sr1);
              mul_r(optr->prop0, (0.5/optr->kappa), optr->prop0, VOLUME / 2);
              mul_r(optr->prop1, (0.5/optr->kappa), optr->prop1, VOLUME / 2);
              convert_eo_to_lexic(lexicfield3, optr->prop0, optr->prop1);
              check1 += scalar_prod_r(lexicfield1,lexicfield3,VOLUME, 1);
              check2 += scalar_prod_i(lexicfield1,lexicfield3,VOLUME, 1);
              diff(lexicfield4,lexicfield3,lexicfield2,VOLUME);
              check4=scalar_prod_r(lexicfield1,lexicfield1,VOLUME, 1);
              check5=square_norm(lexicfield1,VOLUME ,1);

              optr->kappa=kappa_old;
#endif
              if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0)
                {
                  ofs = fopen (filename, "a");
                  ofs_full = fopen (filename_full, "a");
                  fprintf (ofs, "%d %d %d %d   ", traj, operatorid, internum,
                           snum);
                  fprintf (ofs_full, "%d %d %d %d   ", traj, operatorid,
                           internum, snum);
#ifdef CHECK_OPERATOR
                  fprintf( ofs_full, "%e %e %e %e %e %e %e %e %e %e %e %e\n", square1, square2,prodre,prodim,prodreg5,prodimg5,checkg5,check1,check2,check3,check4,check5);
#else
                  /*print all raw data for cross check*/
                  fprintf (ofs_full,
                           "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g   ", k2mu0,
                           k2mu, kappa0, kappa, csw0, csw, rmu0, rmu);
                  fprintf (ofs_full, "%.17g %.17g %.17g %.17g %.17g %.17g\n",
                           square1, square2, prodre, prodim, prodreg5,
                           prodimg5);
                  /*print two and one flavour reweighting log factors*/
                  if (murew)
                    { /* ignoring rounding errors*/
                      fprintf (ofs, "%.17g %.17g %.17g\n",
                               (rmu * rmu - rmu0 * rmu0) * square2,
                               (rmu - rmu0) * prodreg5,
                               (rmu - rmu0) * prodimg5);
                    }
                  if (kapparew)
                    {
                      mdiff = 0.5 / kappa - 0.5 / kappa0;
                      fprintf (ofs, "%.17g %.17g %.17g\n",
                               2.0 * mdiff * prodre + mdiff * mdiff * square2,
                               mdiff * prodre, mdiff * prodim);
                    }
#endif
                  fclose (ofs);
                  fclose (ofs_full);
                }
            }
          else
            { /* end not even odd*/
              set_even_to_zero (optr->sr0);
              set_even_to_zero (optr->prop0);
              random_spinor_field_eo (lexicfield1, reproduce_randomnumber_flag,
                                      RN_GAUSS);
              square1 = square_norm (lexicfield1, VOLUME / 2, 1);

              optr->kappa = kappa0;
              optr->mu = k2mu0;
              update_global_parameters (operatorid);

              optr->applyQm (optr->sr0, lexicfield1);
              assign (optr->prop0, lexicfield1, VOLUME / 2);

              optr->kappa = kappa;
              optr->mu = k2mu;
              /* done in inverter: update_global_parameters(operatorid);*/

              invert_operator_Q (optr->prop0, optr->sr0, operatorid, 0,updateinverter);
              updateinverter=0;

              square2 = square_norm (optr->prop0, VOLUME / 2, 1);

              if (g_mpi_time_rank == 0 && g_proc_coords[0] == 0)
                {
                  ofs = fopen (filename, "a");
                  ofs_full = fopen (filename_full, "a");
                  fprintf (ofs, "%d %d %d %d ", traj, operatorid, internum,
                           snum);
                  fprintf (ofs_full, "%d %d %d %d ", traj, operatorid, internum,
                           snum);
                  fprintf (ofs_full,
                           "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g   ", k2mu0,
                           k2mu, kappa0, kappa, csw0, csw, rmu, rmu0);
                  fprintf (ofs_full, "%.17g %.17g %.17g\n", square1, square2,
                           cswpart);
                  fprintf (ofs, "%.17g\n", square1 - square2 + cswpart);
                  fclose (ofs);
                  fclose (ofs_full);
                }

            }/* end even odd*/

        }/* loop over estimators */

    } /* loop over interpolation steps*/

  etime = gettime ();

  if (g_proc_id == 0 && g_debug_level > 0)
    {
      printf ("REWEIGHTING: measurement done int t/s = %1.4e\n", etime - atime);
    }
  return;
}

void
free_reweighting_parameter (void* par)
{
  reweighting_parameter* param;
  param = (reweighting_parameter*) (par);
  if (param->coeff.el)
    free (param->coeff.el);
  if (param->splitlist.ord)
    free (param->splitlist.ord);
  if (param->splitlist.est)
    free (param->splitlist.est);
  param->coeff.el = NULL;
  param->coeff.s = 0;
  param->splitlist.ord = NULL;
  param->splitlist.est = NULL;
  param->splitlist.s = 0;
}

static void
read_coeff_from_file (vector_list* v)
{
  FILE* file;
  unsigned int l;
  double dat;

  file = fopen ("coeff.dat", "r");
  l = 0;
  dat = 0;

  if (file)
    {
      while (fscanf (file, "%lg ", &dat) > 0)
        {
          l++;
        }
      v->s = l;
      v->el = malloc (l * sizeof(double));
      file = freopen ("coeff.dat", "r", file);
      l = 0;
      while (fscanf (file, "%lg ", &dat) > 0)
        {
          v->el[l++] = dat;
        }
      if (g_debug_level > 3)
        {
          printf (
              "The following coefficients have been read from file coeff.dat:\n");
          for (l = 0; l < v->s; l++)
            {
              printf ("%d %lg\n", l, v->el[l]);
            }
        }
    }
  else
    {
      if (g_proc_id == 0)
        {
          printf ("File coeff.dat not present.\n");
        }
      v->el = NULL;
      v->s = 0;
    }

}

static void
read_splitlist (split_list* list)
{
  FILE* file;
  unsigned int l;
  int dat1, dat2;

  file = fopen ("split.dat", "r");
  l = 0;
  dat1 = dat2 = 0;

  if (file)
    {
      while (fscanf (file, "%d ", &dat1) > 0 && fscanf (file, "%d ", &dat2) > 0)
        {
          l++;
        }
      list->s = l;
      list->ord = malloc (l * sizeof(unsigned int));
      list->est = malloc (l * sizeof(unsigned int));

      file = freopen ("split.dat", "r", file);
      l = 0;
      while (fscanf (file, "%d ", &dat1) > 0 && fscanf (file, "%d ", &dat2) > 0)
        {
          list->ord[l] = dat1;
          list->est[l] = dat2;
          l++;
        }
      if (g_debug_level > 3)
        {
          printf (
              "The following factor splits have been read from file split.dat:\n");
          for (l = 0; l < list->s; l++)
            {
              printf ("%d %d %d\n", l, list->ord[l], list->est[l]);
            }
        }
    }
  else
    {
      if (g_proc_id == 0)
        {
          printf ("File split.dat not present.\n");
        }
      list->ord = NULL;
      list->est = NULL;
      list->s = 0;
    }

}

void
initialize_reweighting_parameter (void** parameter)
{
  reweighting_parameter* param;
  if (!(*parameter))
    {
      (*parameter) = malloc (sizeof(reweighting_parameter));
      param = (reweighting_parameter*) (*parameter);
      param->reweighting_operator = 0;
      param->reweighting_number_sources = 0;
      param->use_evenodd = 0;
      param->k2mu0 = 0.0;
      param->kappa0 = 0.0;
      param->rmu0 = 0.0;
      param->rmu = 0.0;
      param->minev = 1e-7;
      param->maxev = 20.0;
      param->interpolationsteps = 1;
      param->estimatorscheb = 0;
      param->cheborder = 0;
      param->use_cheb = 0;
      param->only_cheb = 0;
      param->coeff.el = NULL;
      param->coeff.s = 0;
      param->splitlist.ord = NULL;
      param->splitlist.est = NULL;
      param->splitlist.s = 0;
      param->evest = 0;
      param->testchebconvergence = 0;
      read_coeff_from_file (&param->coeff);
      read_splitlist (&param->splitlist);
    }
}
