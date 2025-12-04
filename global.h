/***********************************************************************
 *
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Modified by Jenifer Gonzalez Lopez 31.03.2009
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

#ifndef _GLOBAL_H
#define _GLOBAL_H
/***************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 *
 ***************************************************************/
#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#ifdef FIXEDVOLUME
#include "fixed_volume.h"
#endif
#include "su3.h"
#include "su3adj.h"

#include "gettime.h"
#include "misc_types.h"

#define N_CHEBYMAX 49
#define NTILDE_CHEBYMAX 2000

/* size of the extra_masses array for operators using the CGMMS solver */
#define MAX_EXTRA_MASSES 30

#if defined INIT_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

#if ((defined SSE) || (defined SSE2) || (defined SSE3))
#include "sse.h"
#elif defined BGL
#include "bgl.h"
#endif
EXTERN int DUM_DERI, DUM_MATRIX;
EXTERN int NO_OF_SPINORFIELDS;
EXTERN int NO_OF_SPINORFIELDS_32;

EXTERN int DUM_BI_DERI, DUM_BI_SOLVER, DUM_BI_MATRIX;
EXTERN int NO_OF_BISPINORFIELDS;

EXTERN int g_update_gauge_copy;
EXTERN int g_update_gauge_copy_32;
EXTERN int g_relative_precision_flag;
EXTERN int g_strict_residual_check;
EXTERN int g_debug_level;
EXTERN int g_disable_IO_checks;
EXTERN int g_disable_src_IO_checks;

EXTERN tm_mpi_thread_level_t g_mpi_thread_level;

EXTERN tm_timers_t g_timers;

EXTERN int T_global;
#ifndef FIXEDVOLUME
EXTERN int T, L, LX, LY, LZ, VOLUME;
EXTERN int N_PROC_T, N_PROC_X, N_PROC_Y, N_PROC_Z;
EXTERN int RAND, EDGES, VOLUMEPLUSRAND;
EXTERN int TEOSLICE;
EXTERN int SPACEVOLUME, SPACERAND;
#endif

/* translates from lexicographic order to even/odd order */
EXTERN int *g_lexic2eo;
/* translates from even/odd order to lexicograhic order  */
EXTERN int *g_eo2lexic;
EXTERN int *g_lexic2eosub;
EXTERN int g_sloppy_precision_flag;
EXTERN int g_sloppy_precision;

EXTERN int ****g_ipt;
EXTERN int **g_iup;
EXTERN int **g_idn;
EXTERN int **g_iup_eo; /* NEW GIUPDNEO */
EXTERN int **g_idn_eo;
EXTERN int **g_coord;
EXTERN int *g_hi;

EXTERN int *g_field_z_ipt_even;
EXTERN int *g_field_z_ipt_odd;

EXTERN spinor **g_spinor_field;
EXTERN spinor32 **g_spinor_field32;

EXTERN bispinor **g_bispinor_field;
EXTERN spinor *g_tbuff;

/* Index independent geometry */

EXTERN int *g_field_z_ipt_even;
EXTERN int *g_field_z_ipt_odd;
EXTERN int *g_field_z_disp_even_dn;
EXTERN int *g_field_z_disp_even_up;
EXTERN int *g_field_z_disp_odd_dn;
EXTERN int *g_field_z_disp_odd_up;


/* IF PHMC  */
EXTERN spinor **g_chi_up_spinor_field;
EXTERN spinor **g_chi_dn_spinor_field;
EXTERN int g_running_phmc;
/* End IF PHMC  */

EXTERN su3 **g_gauge_field;
EXTERN su3_32 **g_gauge_field_32;
#ifdef _USE_HALFSPINOR
EXTERN su3 ***g_gauge_field_copy;
EXTERN su3_32 ***g_gauge_field_copy_32;
#elif (defined _USE_TSPLITPAR)
EXTERN su3 **g_gauge_field_copyt;
EXTERN su3 **g_gauge_field_copys;
#else
EXTERN su3 **g_gauge_field_copy;
EXTERN su3_32 **g_gauge_field_copy_32;
#endif

EXTERN su3adj **moment;
EXTERN su3adj **df0;
EXTERN su3adj **ddummy, **debug_derivative;

EXTERN int count00, count01, count10, count11, count20, count21;
EXTERN double g_kappa, g_c_sw, g_beta;
EXTERN double g_mu, g_mu1, g_mu2, g_mu3, g_shift;
EXTERN double g_rgi_C0, g_rgi_C1;

/* Parameters for non-degenrate case */
EXTERN double g_mubar, g_epsbar;
EXTERN int g_use_clover_flag;

/* MPI information */
EXTERN int g_proc_id, g_nproc, g_stdio_proc, g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z, g_cart_id;
EXTERN int g_proc_coords[4];
EXTERN int g_dbw2rand;
EXTERN int g_mpi_time_rank;
EXTERN int g_mpi_SV_rank;
EXTERN int g_mpi_z_rank;
EXTERN int g_mpi_ST_rank;
EXTERN int g_nb_list[8];

/* Variables for exposu3 */
EXTERN int g_exposu3_no_c;
EXTERN double *g_exposu3_c;

/* OpenMP Kahan accumulation arrays */
EXTERN _Complex double *g_omp_acc_cp;
EXTERN double *g_omp_acc_re;

/* Deflation information */
EXTERN int g_dflgcr_flag;
EXTERN int g_N_s;
EXTERN int *index_block_eo;
EXTERN int Msap_precon;
EXTERN int NiterMsap;
EXTERN int NcycleMsap;
EXTERN int NiterMsap_dflgen;
EXTERN int NcycleMsap_dflgen;
EXTERN int NsmoothMsap_dflgen;
EXTERN int usePL;
EXTERN int little_solver;
EXTERN int little_evenodd;
EXTERN int little_gmres_m_parameter;
EXTERN double little_solver_low_prec;
EXTERN double little_solver_high_prec;
EXTERN int little_solver_max_iter;

#ifdef TM_USE_MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1, req2, req3, req4;
EXTERN MPI_Comm g_cart_grid;
EXTERN MPI_Comm g_mpi_time_slices;
EXTERN MPI_Comm g_mpi_SV_slices;
EXTERN MPI_Comm g_mpi_z_slices;
EXTERN MPI_Comm g_mpi_ST_slices;

/* the next neighbours for MPI */
EXTERN int g_nb_x_up, g_nb_x_dn;
EXTERN int g_nb_y_up, g_nb_y_dn;
EXTERN int g_nb_t_up, g_nb_t_dn;
EXTERN int g_nb_z_up, g_nb_z_dn;

#endif

#ifdef TM_USE_OMP
EXTERN int omp_num_threads;
#endif

/* something to evaluate time elaps */
EXTERN double DeltaTtot, DeltaTcd, DeltaTev;
EXTERN int counter_Spsi;
/* end of the something ... */

EXTERN void *g_precWS;

#ifdef WITHLAPH
/* Jacobi operator per Laplacian Heaviside (LapH) */
EXTERN su3_vector **g_jacobi_field;
EXTERN int gI_0_0_0, gI_L_0_0, gI_Lm1_0_0, gI_m1_0_0, gI_0_L_0, gI_0_Lm1_0, gI_0_m1_0, gI_0_0_L,
    gI_0_0_Lm1, gI_0_0_m1;
EXTERN int tempT, tempV, tempR;
EXTERN int **g_iup3d;
EXTERN int **g_idn3d;
#endif

/* keeping track of what the gauge, clover and inverse clover
 * field contain in order to avoid unnecessary inversions
 * of the latter */
EXTERN tm_GaugeState_t g_gauge_state;
EXTERN tm_GaugeState_t g_gauge_state_32;
EXTERN tm_CloverState_t g_clover_state;
EXTERN tm_CloverState_t g_clover_state_32;
EXTERN tm_CloverInverseState_t g_clover_inverse_state;
EXTERN tm_CloverInverseState_t g_clover_inverse_state_32;

#undef EXTERN
/* #undef ALIGN */

void fatal_error(char const *error, char const *function);

#endif

/*
 * Comments: generic macro for swapping values or pointers.
 * We use memcpy because is optimal when the amount to copy is known at compilation time.
 * "sizeof(x) == sizeof(y) ? (signed)sizeof(x) : -1" is a compile time check that the types are
 * compatible.
 */
#define SWAP(x, y)                                                            \
  do {                                                                        \
    unsigned char swap_temp[sizeof(x) == sizeof(y) ? (signed)sizeof(x) : -1]; \
    memcpy(swap_temp, &y, sizeof(x));                                         \
    memcpy(&y, &x, sizeof(x));                                                \
    memcpy(&x, swap_temp, sizeof(x));                                         \
  } while (0)
