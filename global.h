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
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#ifdef MPI
#  include <mpi.h>
#endif
#ifdef FIXEDVOLUME
#  include "fixed_volume.h"
#endif
#include "su3.h"
#include "su3adj.h"
//#  include <tormpi_export.h>

#define N_CHEBYMAX 49
#define NTILDE_CHEBYMAX 2000

#if defined MAIN_PROGRAM
#  define EXTERN
#else
#  define EXTERN extern
#endif

#if ((defined SSE)||(defined SSE2)||(defined SSE3))
#  include "sse.h"
#elif defined BGL
# include "bgl.h"
#else
#  define ALIGN
#endif

EXTERN int DUM_DERI, DUM_SOLVER, DUM_MATRIX;
EXTERN int NO_OF_SPINORFIELDS;

EXTERN int DUM_BI_DERI, DUM_BI_SOLVER, DUM_BI_MATRIX;
EXTERN int NO_OF_BISPINORFIELDS;

EXTERN int g_update_gauge_copy;
EXTERN int g_update_gauge_energy;
EXTERN int g_update_rectangle_energy;
EXTERN int g_relative_precision_flag;
EXTERN int g_debug_level;
EXTERN int g_disable_IO_checks;

EXTERN int T_global;
#ifndef FIXEDVOLUME
EXTERN int T, L, LX, LY, LZ, VOLUME;
EXTERN int N_PROC_T, N_PROC_X, N_PROC_Y, N_PROC_Z;
EXTERN int RAND, EDGES, VOLUMEPLUSRAND;
EXTERN int TEOSLICE;
EXTERN int SPACEVOLUME, SPACERAND;
#endif

/* translates from lexicographic order to even/odd order */
EXTERN int * g_lexic2eo;
/* translates from even/odd order to lexicograhic order  */
EXTERN int * g_eo2lexic;
EXTERN int * g_lexic2eosub;
EXTERN int g_sloppy_precision_flag;
EXTERN int g_sloppy_precision;

EXTERN int **** g_ipt;
EXTERN int ** g_iup;
EXTERN int ** g_idn;
EXTERN int ** g_iup_eo; /* NEW GIUPDNEO */
EXTERN int ** g_idn_eo;
EXTERN int ** g_coord;


EXTERN int * g_field_z_ipt_even;
EXTERN int * g_field_z_ipt_odd;

EXTERN spinor ** g_spinor_field;

EXTERN bispinor ** g_bispinor_field;
EXTERN spinor * g_tbuff;

/* Index independent geometry */

EXTERN int * g_field_z_ipt_even;
EXTERN int * g_field_z_ipt_odd;
EXTERN int * g_field_z_disp_even_dn;
EXTERN int * g_field_z_disp_even_up;
EXTERN int * g_field_z_disp_odd_dn;
EXTERN int * g_field_z_disp_odd_up;

/* this if statement will be removed in future and _INDEX_INDEP_GEOM will be the default */
#ifdef _INDEX_INDEP_GEOM
EXTERN int g_1st_t_int_dn,g_1st_t_int_up,g_1st_t_ext_dn,g_1st_t_ext_up;
EXTERN int g_1st_x_int_dn,g_1st_x_int_up,g_1st_x_ext_dn,g_1st_x_ext_up;
EXTERN int g_1st_y_int_dn,g_1st_y_int_up,g_1st_y_ext_dn,g_1st_y_ext_up;
EXTERN int g_1st_z_int_dn,g_1st_z_int_up,g_1st_z_ext_dn,g_1st_z_ext_up;
EXTERN int gI_0_0_0_0,gI_L_0_0_0,gI_Lm1_0_0_0,gI_m1_0_0_0,gI_p1_0_0_0,gI_Lp1_0_0_0,gI_Lm2_0_0_0,gI_m2_0_0_0;
EXTERN int gI_0_L_0_0,gI_0_Lm1_0_0,gI_0_m1_0_0,gI_0_p1_0_0,gI_0_Lp1_0_0,gI_0_Lm2_0_0,gI_0_m2_0_0,gI_L_L_0_0;
EXTERN int gI_Lm1_L_0_0,gI_m1_L_0_0,gI_p1_L_0_0,gI_Lp1_L_0_0,gI_Lm2_L_0_0,gI_m2_L_0_0,gI_L_Lp1_0_0,gI_Lm1_Lp1_0_0;
EXTERN int gI_m1_Lp1_0_0,gI_0_0_L_0,gI_0_0_Lm1_0,gI_0_0_m1_0,gI_0_0_p1_0,gI_0_0_Lp1_0,gI_0_0_Lm2_0,gI_0_0_m2_0;
EXTERN int gI_0_L_L_0,gI_0_Lm1_L_0,gI_0_m1_L_0,gI_L_0_L_0,gI_L_0_Lm1_0,gI_L_0_m1_0,gI_0_p1_L_0,gI_0_Lp1_L_0;
EXTERN int gI_0_Lm2_L_0,gI_0_m2_L_0,gI_0_L_Lp1_0,gI_0_Lm1_Lp1_0,gI_0_m1_Lp1_0,gI_Lp1_0_L_0,gI_Lp1_0_Lm1_0;
EXTERN int gI_Lp1_0_m1_0,gI_L_0_p1_0,gI_L_0_Lp1_0,gI_L_0_Lm2_0,gI_L_0_m2_0,gI_0_0_0_L,gI_0_0_0_Lm1,gI_0_0_0_m1;
EXTERN int gI_0_0_0_p1,gI_0_0_0_Lp1,gI_0_0_0_Lm2,gI_0_0_0_m2,gI_0_L_0_L,gI_0_Lm1_0_L,gI_0_m1_0_L,gI_L_0_0_L;
EXTERN int gI_L_0_0_Lm1,gI_L_0_0_m1,gI_0_L_0_L,gI_0_Lm1_0_L,gI_0_m1_0_L,gI_Lp1_0_0_L,gI_Lp1_0_0_Lm1,gI_Lp1_0_0_m1;
EXTERN int gI_L_0_0_p1,gI_L_0_0_Lp1,gI_L_0_0_Lm2,gI_L_0_0_m2,gI_0_L_0_Lp1,gI_0_Lm1_0_Lp1,gI_0_m1_0_Lp1,gI_0_p1_0_L;
EXTERN int gI_0_Lp1_0_L,gI_0_Lm2_0_L,gI_0_m2_0_L,gI_0_0_L_Lp1,gI_0_0_Lm1_Lp1,gI_0_0_m1_Lp1,gI_0_0_p1_L;
EXTERN int gI_0_0_Lp1_L,gI_0_0_Lm2_L,gI_0_0_m2_L,gI_Lp1_m1_0_0,gI_m2_m1_0_0,gI_m2_0_L_0,gI_m2_0_m1_0,gI_0_Lp1_m1_0;
EXTERN int gI_0_m2_m1_0,gI_m2_0_0_L,gI_m2_0_0_m1,gI_0_Lp1_0_m1,gI_0_m2_0_m1,gI_0_0_Lp1_m1,gI_0_0_m2_m1,gI_m1_0_0_m2;
EXTERN int gI_0_0_L_L, gI_0_0_m1_L, gI_0_0_Lm1_L;

# ifdef _USE_HALFSPINOR
EXTERN int g_HS_shift_t,g_HS_shift_x,g_HS_shift_y,g_HS_shift_z;
# endif

# ifdef _USE_TSPLITPAR
EXTERN int ** g_field_zt_disp_even_dn;
EXTERN int ** g_field_zt_disp_even_up;
EXTERN int ** g_field_zt_disp_odd_dn;
EXTERN int ** g_field_zt_disp_odd_up;
EXTERN int ** g_1st_eot;
EXTERN int * g_1st_xt_int_dn;
EXTERN int * g_1st_xt_int_up;
EXTERN int * g_1st_xt_ext_dn;
EXTERN int * g_1st_xt_ext_up;
EXTERN int * g_1st_yt_int_dn;
EXTERN int * g_1st_yt_int_up;
EXTERN int * g_1st_yt_ext_dn;
EXTERN int * g_1st_yt_ext_up;
EXTERN int * g_1st_zt_int_dn;
EXTERN int * g_1st_zt_int_up;
EXTERN int * g_1st_zt_ext_dn;
EXTERN int * g_1st_zt_ext_up;
# endif
#endif /* _INDEX_INDEP_GEOM */ 

/* IF PHMC  */
EXTERN spinor ** g_chi_up_spinor_field;
EXTERN spinor ** g_chi_dn_spinor_field;
EXTERN int g_running_phmc;
/* End IF PHMC  */

EXTERN su3 ** g_gauge_field;
#ifdef _USE_HALFSPINOR
EXTERN su3 *** g_gauge_field_copy;
#elif (defined _USE_TSPLITPAR )
EXTERN su3 ** g_gauge_field_copyt;
EXTERN su3 ** g_gauge_field_copys;
#else
EXTERN su3 ** g_gauge_field_copy;
#endif

/*for temporalgauge in GPU part*/
EXTERN su3 ** g_tempgauge_field;

EXTERN su3adj ** moment;
EXTERN su3adj ** df0;
EXTERN su3adj ** ddummy;

EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double g_kappa, g_c_sw, g_ka_csw_8, g_beta;
EXTERN double g_mu, g_mu1, g_mu2, g_mu3;
EXTERN double g_rgi_C0, g_rgi_C1;

/*************************/
/* SF definitions */
EXTERN double g_eta;                          /* background field parameter */
EXTERN double g_Ct, g_Cs;                     /* plaquette part */
EXTERN double g_C1ss, g_C1tss, g_C1tts;       /* rectangle part */
EXTERN int g_Tbsf;                            /* it sets at which time slice I want to put the SF b.c. (end point)
                                                 T = lattice time extent */
/* variables specifying the value of t,x,y,z for each lattice site ix */
EXTERN int* g_t;
EXTERN int* g_x;
EXTERN int* g_y;
EXTERN int* g_z;
EXTERN int g_sf_inc_wrap_sq;
/* end of SF definitions */
/*************************/

/* Parameters for non-degenrate case */
EXTERN double g_acc_Ptilde;
EXTERN double g_acc_Hfin;
EXTERN int g_rec_ev;
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

/* Deflation information */
EXTERN int g_dflgcr_flag;
EXTERN int g_N_s;
EXTERN int * index_block_eo;

#ifdef MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
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

/* something to evaluate time elaps */
EXTERN double DeltaTtot, DeltaTcd, DeltaTev;
EXTERN int counter_Spsi;
/* end of the something ... */

/* Info for CGMMS solver */
EXTERN int g_no_extra_masses;
EXTERN double g_extra_masses[30];

EXTERN int ITER_MAX_BCG;
EXTERN int ITER_MAX_CG;

EXTERN void* g_precWS;

#ifdef WITHLAPH
/* Jacobi operator per Laplacian Heaviside (LapH) */
EXTERN su3_vector ** g_jacobi_field;
EXTERN int gI_0_0_0, gI_L_0_0, gI_Lm1_0_0, gI_m1_0_0, gI_0_L_0, gI_0_Lm1_0, gI_0_m1_0, gI_0_0_L, gI_0_0_Lm1, gI_0_0_m1;
EXTERN int tempT,tempV,tempR;
EXTERN int ** g_iup3d;
EXTERN int ** g_idn3d;
#endif

#undef EXTERN
/* #undef ALIGN */

void fatal_error(char const *error, char const *function);

#ifdef MAIN_PROGRAM
#endif

#endif

