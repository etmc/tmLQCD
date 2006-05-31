/* $Id$ */
#ifndef _GLOBAL_H
#define _GLOBAL_H
/***************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * Author: Martin Luescher <luscher@mail.desy.de>
 * Date: 16.03.2001
 *
 * Adapted for the HMC-Program by M. Hasenbusch 2002
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
#include"su3.h"
#include"su3adj.h"

#define DUM_DERI 6
#define DUM_SOLVER (DUM_DERI+7)
#define DUM_MATRIX (DUM_SOLVER+6)
/* if you want to include bicgstabell */
/* #define DUM_MATRIX (DUM_SOLVER+11) */
#define NO_OF_SPINORFIELDS (DUM_MATRIX+6)

/* Here you can define antiperiodic  */
/* boundary conditions with e.g.     */
/* #define X1 1.  (in x-direction)          */
#define X1 0.
#define X2 0.
#define X3 0.
#define EPS_SQ0  1.0e-16
#define EPS_SQ1  1.0e-7
#define EPS_SQ2  1.0e-5
#define EPS_SQ3  1.0e-3
#define tiny_t  1.0e-20

#define N_CHEBYMAX 2000


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

EXTERN double g_eps_sq_force, g_eps_sq_acc;
EXTERN double g_eps_sq_force1, g_eps_sq_force2, g_eps_sq_force3;
EXTERN double g_eps_sq_acc1, g_eps_sq_acc2, g_eps_sq_acc3;
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
/* translates from even/odd orderto lexicograhic order  */    
EXTERN int * g_eo2lexic;
EXTERN int * g_lexic2eosub;

EXTERN int **** g_ipt;
EXTERN int ** g_iup;
EXTERN int ** g_idn;
EXTERN int ** g_halfpt;

EXTERN int * g_field_z_ipt_even;
EXTERN int * g_field_z_ipt_odd;

EXTERN spinor ** g_spinor_field;

EXTERN bispinor ** g_bispinor_field;

EXTERN spinor ** g_csg_field[4];
EXTERN int * g_csg_index_array[4];
EXTERN int g_csg_N[8];

EXTERN su3 ** g_gauge_field;
EXTERN su3 ** g_gauge_field_copy;

EXTERN su3adj ** moment;
EXTERN su3adj ** df0;
EXTERN su3adj ** ddummy;

EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double g_kappa, g_c_sw, g_ka_csw_8, g_beta;
EXTERN double g_rgi_C0, g_rgi_C1;
EXTERN double g_mu, g_mu1, g_mu2, g_mu3;
/* Parameters for non-degenrate case */
EXTERN double g_mubar, g_epsbar;
EXTERN int g_use_clover_flag, g_nr_of_psf;

/* MPI information */
EXTERN int g_proc_id, g_nproc, g_stdio_proc, g_nproc_t, g_nproc_x, g_nproc_y, g_nproc_z, g_cart_id;
EXTERN int g_proc_coords[4];
EXTERN int g_dbw2rand;
#ifdef MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
EXTERN MPI_Comm g_cart_grid;

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

EXTERN int ITER_MAX_BCG;
EXTERN int ITER_MAX_CG;

#undef EXTERN
/* #undef ALIGN */

#ifdef MAIN_PROGRAM
static char const g_rcsid[] = "$Id$";
#endif

#endif

