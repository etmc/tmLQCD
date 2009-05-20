/***********************************************************************
 * $Id$
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

EXTERN int T_global;
#ifndef FIXEDVOLUME
EXTERN int T, L, LX, LY, LZ, VOLUME;
EXTERN int N_PROC_T, N_PROC_X, N_PROC_Y, N_PROC_Z;
EXTERN int RAND, EDGES, VOLUMEPLUSRAND;
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

EXTERN int * g_field_z_ipt_even;
EXTERN int * g_field_z_ipt_odd;

EXTERN spinor ** g_spinor_field;

EXTERN bispinor ** g_bispinor_field;

/* IF PHMC  */
EXTERN spinor ** g_chi_up_spinor_field;
EXTERN spinor ** g_chi_dn_spinor_field;
EXTERN int g_running_phmc;
/* End IF PHMC  */

EXTERN su3 ** g_gauge_field;
#ifdef _USE_HALFSPINOR
EXTERN su3 *** g_gauge_field_copy;
#else
EXTERN su3 ** g_gauge_field_copy;
#endif

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
                                                 T = lattice time extent set by Carsten */
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
EXTERN int g_nb_list[8];

/* Deflation information */
EXTERN int g_dflgcr_flag;
EXTERN int g_N_s;

#ifdef MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
EXTERN MPI_Comm g_cart_grid;
EXTERN MPI_Comm g_mpi_time_slices;
EXTERN MPI_Comm g_mpi_SV_slices;

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

#undef EXTERN
/* #undef ALIGN */

#ifdef MAIN_PROGRAM
static char const g_rcsid[] = "$Id$";
#endif

#endif

