/***********************************************************************
 *
 * Copyright (C) 1008 Carsten Urbach
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
#include "boundary.h"
#include "deriv_Sb.h"
#include "deriv_Sb_D_psi.h"
#include "detratio_monomial.h"
#include "gamma.h"
#include "gettime.h"
#include "global.h"
#include "linalg_eo.h"
#include "misc_types.h"
#include "monomial/monomial.h"
#include "operator/Hopping_Matrix.h"
#include "operator/tm_operators.h"
#include "read_input.h"
#include "solver/chrono_guess.h"
#include "solver/monomial_solve.h"
#include "solver/solver.h"
#include "start.h"
#include "su3.h"

/* think about chronological solver ! */

void detratio_derivative(const int no, hamiltonian_field_t* const hf) {
  monomial* mnl = &monomial_list[no];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);

  mnl->forcefactor = 1.;

  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);
  g_mu = mnl->mu2;
  g_kappa = mnl->kappa2;
  boundary(g_kappa);

  if (mnl->even_odd_flag) {
    /*
     *  this is being run in case there is even/odd preconditioning
     */
    /*********************************************************************
     *
     * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
     * mu1 and mu2 are set according to the monomial
     *
     *********************************************************************/
    /* First term coming from the second field */
    /* Multiply with W_+ */
    tm_stopwatch_push(&g_timers, "Qtm_plus_psi", "");
    Qtm_plus_psi(mnl->w_fields[2], mnl->pf);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu;
    g_kappa = mnl->kappa;
    boundary(g_kappa);

    /* Invert Q_{+} Q_{-} */
    /* X_W -> w_fields[1] */
    chrono_guess(mnl->w_fields[1], mnl->w_fields[2], mnl->csg_field, mnl->csg_index_array,
                 mnl->csg_N, mnl->csg_n, VOLUME / 2, &Qtm_pm_psi);

    mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->w_fields[2], mnl->solver_params,
                                   mnl->maxiter, mnl->forceprec, g_relative_precision_flag,
                                   VOLUME / 2, &Qtm_pm_psi, mnl->solver);

    chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array, mnl->csg_N,
                        &mnl->csg_n, VOLUME / 2);
    /* Y_W -> w_fields[0]  */
    tm_stopwatch_push(&g_timers, "Qtm_minus_psi", "");
    Qtm_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    tm_stopwatch_push(&g_timers, "H_eo_tm_inv_psi", "");
    H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1.);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor);

    /* to get the even sites of Y */
    tm_stopwatch_push(&g_timers, "H_eo_tm_inv_psi", "");
    H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf, mnl->forcefactor);

    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(mnl->w_fields[0], -1., mnl->pf, VOLUME / 2);

    g_mu = mnl->mu2;
    g_kappa = mnl->kappa2;
    boundary(g_kappa);

    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    tm_stopwatch_push(&g_timers, "H_eo_tm_inv_psi", "");
    H_eo_tm_inv_psi(mnl->w_fields[2], mnl->w_fields[1], EO, -1.);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, mnl->w_fields[0], mnl->w_fields[2], hf, mnl->forcefactor);

    /* to get the even sites of Y */
    tm_stopwatch_push(&g_timers, "H_eo_tm_inv_psi", "");
    H_eo_tm_inv_psi(mnl->w_fields[3], mnl->w_fields[0], EO, +1);
    tm_stopwatch_pop(&g_timers, 0, 1, "");
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, mnl->w_fields[3], mnl->w_fields[1], hf, mnl->forcefactor);
  } else { /* no even/odd preconditioning */
    /*********************************************************************
     *
     * This term is det((Q^2 + \mu_1^2)/(Q^2 + \mu_2^2))
     * mu1 and mu2 are set according to the monomial
     *
     *********************************************************************/
    /* First term coming from the second field */
    /* Multiply with W_+ */
    tm_stopwatch_push(&g_timers, "Q_plus_psi", "");
    Q_plus_psi(mnl->w_fields[2], mnl->pf);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu;
    g_kappa = mnl->kappa;
    boundary(g_kappa);

    if ((mnl->solver == CG) || (mnl->solver == MIXEDCG) || (mnl->solver == RGMIXEDCG) ||
        (mnl->solver == MG)) {
      /* If CG is used anyhow */
      /*       gamma5(mnl->w_fields[1], mnl->w_fields[2], VOLUME/2); */
      /* Invert Q_{+} Q_{-} */
      /* X_W -> w_fields[1] */
      chrono_guess(mnl->w_fields[1], mnl->w_fields[2], mnl->csg_field, mnl->csg_index_array,
                   mnl->csg_N, mnl->csg_n, VOLUME / 2, &Q_pm_psi);
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->w_fields[2], mnl->solver_params,
                                     mnl->maxiter, mnl->forceprec, g_relative_precision_flag,
                                     VOLUME, &Q_pm_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array, mnl->csg_N,
                          &mnl->csg_n, VOLUME / 2);

      /* Y_W -> w_fields[0]  */
      tm_stopwatch_push(&g_timers, "Q_plus_psi", "");
      Q_minus_psi(mnl->w_fields[0], mnl->w_fields[1]);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    } else {
      /* Invert first Q_+ */
      /* Y_o -> w_fields[0]  */

      chrono_guess(mnl->w_fields[0], mnl->w_fields[2], mnl->csg_field, mnl->csg_index_array,
                   mnl->csg_N, mnl->csg_n, VOLUME / 2, &Q_plus_psi);
      gamma5(mnl->w_fields[0], mnl->w_fields[0], VOLUME);
      mnl->iter1 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[2], mnl->solver_params,
                                     mnl->maxiter, mnl->forceprec, g_relative_precision_flag,
                                     VOLUME, Q_plus_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[0], mnl->csg_field, mnl->csg_index_array, mnl->csg_N,
                          &mnl->csg_n, VOLUME / 2);

      /* Now Q_- */
      /* X_o -> w_fields[1] */

      chrono_guess(mnl->w_fields[1], mnl->w_fields[0], mnl->csg_field2, mnl->csg_index_array2,
                   mnl->csg_N2, mnl->csg_n2, VOLUME / 2, &Q_minus_psi);
      gamma5(mnl->w_fields[1], mnl->w_fields[1], VOLUME);
      mnl->iter1 += solve_degenerate(mnl->w_fields[1], mnl->w_fields[0], mnl->solver_params,
                                     mnl->maxiter, mnl->forceprec, g_relative_precision_flag,
                                     VOLUME, Q_minus_psi, mnl->solver);
      chrono_add_solution(mnl->w_fields[1], mnl->csg_field2, mnl->csg_index_array2, mnl->csg_N2,
                          &mnl->csg_n2, VOLUME / 2);
    }

    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(mnl->w_fields[0], mnl->w_fields[1], hf, mnl->forcefactor);

    g_mu = mnl->mu2;
    g_kappa = mnl->kappa2;
    boundary(g_kappa);

    /* Second term coming from the second field */
    /* The sign is opposite!! */
    mul_r(mnl->w_fields[0], -1., mnl->pf, VOLUME);

    /* \delta Q sandwitched by Y^\dagger and X */
    deriv_Sb_D_psi(mnl->w_fields[0], mnl->w_fields[1], hf, mnl->forcefactor);
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}

void detratio_heatbath(const int id, hamiltonian_field_t* const hf) {
  monomial* mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);

  g_mu = mnl->mu;
  g_kappa = mnl->kappa;
  boundary(g_kappa);

  mnl->csg_n = 0;
  mnl->csg_n2 = 0;
  mnl->iter0 = 0;
  mnl->iter1 = 0;
  if (mnl->even_odd_flag) {
    tm_stopwatch_push(&g_timers, "random_energy0", "");
    random_spinor_field_eo(mnl->w_fields[0], mnl->rngrepro, RN_GAUSS);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME / 2, 1);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    tm_stopwatch_push(&g_timers, "Qp", "");
    mnl->Qp(mnl->w_fields[1], mnl->w_fields[0]);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu2;
    g_kappa = mnl->kappa2;
    boundary(g_kappa);

    zero_spinor_field(mnl->w_fields[0], VOLUME / 2);
    if (mnl->HB_solver == MG || mnl->HB_solver == BICGSTAB) {
      mnl->iter0 = solve_degenerate(mnl->pf, mnl->w_fields[1], mnl->HB_solver_params,
                                    mnl->HB_maxiter, mnl->accprec, g_relative_precision_flag,
                                    VOLUME / 2, mnl->Qp, mnl->HB_solver);
      chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array, mnl->csg_N, &mnl->csg_n,
                          VOLUME / 2);
    } else {
      mnl->iter0 = solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->HB_solver_params,
                                    mnl->HB_maxiter, mnl->accprec, g_relative_precision_flag,
                                    VOLUME / 2, mnl->Qsq, mnl->HB_solver);
      tm_stopwatch_push(&g_timers, "Qm", "");
      mnl->Qm(mnl->pf, mnl->w_fields[0]);
      tm_stopwatch_pop(&g_timers, 0, 1, "");

      chrono_add_solution(mnl->w_fields[0], mnl->csg_field, mnl->csg_index_array, mnl->csg_N,
                          &mnl->csg_n, VOLUME / 2);
    }
  } else {
    tm_stopwatch_push(&g_timers, "random_energy0", "");
    random_spinor_field_lexic(mnl->w_fields[0], mnl->rngrepro, RN_GAUSS);
    mnl->energy0 = square_norm(mnl->w_fields[0], VOLUME, 1);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    tm_stopwatch_push(&g_timers, "Q_plus_psi", "");
    Q_plus_psi(mnl->w_fields[1], mnl->w_fields[0]);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu2;
    g_kappa = mnl->kappa2;
    boundary(g_kappa);

    zero_spinor_field(mnl->pf, VOLUME);
    if ((mnl->HB_solver == CG) || (mnl->HB_solver == MIXEDCG) || (mnl->HB_solver == RGMIXEDCG)) {
      mnl->iter0 = solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->HB_solver_params,
                                    mnl->HB_maxiter, mnl->accprec, g_relative_precision_flag,
                                    VOLUME, Q_pm_psi, mnl->HB_solver);
      tm_stopwatch_push(&g_timers, "Q_minus_psi", "");
      Q_minus_psi(mnl->pf, mnl->w_fields[0]);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
      chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array, mnl->csg_N, &mnl->csg_n,
                          VOLUME / 2);
      // FIXME why only MG here?
    } else if (mnl->HB_solver == MG) {
      mnl->iter0 = solve_degenerate(mnl->pf, mnl->w_fields[1], mnl->HB_solver_params,
                                    mnl->HB_maxiter, mnl->accprec, g_relative_precision_flag,
                                    VOLUME, Q_plus_psi, mnl->HB_solver);
      chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array, mnl->csg_N, &mnl->csg_n,
                          VOLUME / 2);
    } else {
      mnl->iter0 = solve_degenerate(mnl->pf, mnl->w_fields[1], mnl->HB_solver_params,
                                    mnl->HB_maxiter, mnl->accprec, g_relative_precision_flag,
                                    VOLUME, Q_plus_psi, mnl->HB_solver);
      chrono_add_solution(mnl->pf, mnl->csg_field, mnl->csg_index_array, mnl->csg_N, &mnl->csg_n,
                          VOLUME / 2);
      chrono_add_solution(mnl->pf, mnl->csg_field2, mnl->csg_index_array2, mnl->csg_N2,
                          &mnl->csg_n2, VOLUME / 2);
    }
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  if (g_proc_id == 0) {
    if (g_debug_level > 3) {
      printf("called detratio_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return;
}

double detratio_acc(const int id, hamiltonian_field_t* const hf) {
  monomial* mnl = &monomial_list[id];
  tm_stopwatch_push(&g_timers, __func__, mnl->name);
  int save_sloppy = g_sloppy_precision_flag;

  mnl_backup_restore_globals(TM_BACKUP_GLOBALS);

  g_mu = mnl->mu2;
  g_kappa = mnl->kappa2;
  boundary(g_kappa);

  if (even_odd_flag) {
    tm_stopwatch_push(&g_timers, "Qp", "");
    mnl->Qp(mnl->w_fields[1], mnl->pf);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu;
    g_kappa = mnl->kappa;
    boundary(g_kappa);

    g_sloppy_precision_flag = 0;
    if (mnl->solver == MG || mnl->solver == BICGSTAB) {
      chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
                   mnl->csg_N, mnl->csg_n, VOLUME / 2, mnl->Qp);
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params,
                                     mnl->maxiter, mnl->accprec, g_relative_precision_flag,
                                     VOLUME / 2, mnl->Qp, mnl->solver);
      /* Compute the energy contr. from second field */
      tm_stopwatch_push(&g_timers, "energy1", "");
      mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME / 2, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    } else {
      chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
                   mnl->csg_N, mnl->csg_n, VOLUME / 2, mnl->Qsq);
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params,
                                     mnl->maxiter, mnl->accprec, g_relative_precision_flag,
                                     VOLUME / 2, mnl->Qsq, mnl->solver);
      tm_stopwatch_push(&g_timers, "Qm", "");
      mnl->Qm(mnl->w_fields[1], mnl->w_fields[0]);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
      /* Compute the energy contr. from second field */
      tm_stopwatch_push(&g_timers, "energy1", "");
      mnl->energy1 = square_norm(mnl->w_fields[1], VOLUME / 2, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    }
    g_sloppy_precision_flag = save_sloppy;
  } else {
    tm_stopwatch_push(&g_timers, "Q_plus_psi", "");
    Q_plus_psi(mnl->w_fields[1], mnl->pf);
    tm_stopwatch_pop(&g_timers, 0, 1, "");

    g_mu = mnl->mu;
    g_kappa = mnl->kappa;
    boundary(g_kappa);

    chrono_guess(mnl->w_fields[0], mnl->w_fields[1], mnl->csg_field, mnl->csg_index_array,
                 mnl->csg_N, mnl->csg_n, VOLUME / 2, &Q_plus_psi);
    g_sloppy_precision_flag = 0;
    if ((mnl->solver == CG) || (mnl->solver == MIXEDCG)) {
      mnl->iter0 +=
          solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params, mnl->maxiter,
                           mnl->accprec, g_relative_precision_flag, VOLUME, &Q_pm_psi, mnl->solver);
      tm_stopwatch_push(&g_timers, "Q_minus_psi", "");
      Q_minus_psi(mnl->w_fields[1], mnl->w_fields[0]);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
      /* Compute the energy contr. from second field */
      tm_stopwatch_push(&g_timers, "energy1", "");
      mnl->energy1 = square_norm(mnl->w_fields[1], VOLUME, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    } else {
      mnl->iter0 += solve_degenerate(mnl->w_fields[0], mnl->w_fields[1], mnl->solver_params,
                                     mnl->maxiter, mnl->accprec, g_relative_precision_flag, VOLUME,
                                     Q_plus_psi, mnl->solver);

      /* Compute the energy contr. from second field */
      tm_stopwatch_push(&g_timers, "energy1", "");
      mnl->energy1 = square_norm(mnl->w_fields[0], VOLUME, 1);
      tm_stopwatch_pop(&g_timers, 0, 1, "");
    }
    g_sloppy_precision_flag = save_sloppy;
  }
  mnl_backup_restore_globals(TM_RESTORE_GLOBALS);
  if (g_proc_id == 0) {
    if (g_debug_level > 3) {
      printf("called detratio_acc for id %d dH = %1.10e\n", id, mnl->energy1 - mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, "");
  return (mnl->energy1 - mnl->energy0);
}
