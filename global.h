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
#include <mpi.h>
#endif
#include"su3.h"
#include"su3adj.h"

#define T   2
#define L   4

#ifndef PARALLELXT
#define N_PROC_X 1
#else
#define N_PROC_X 2
#endif

#define LX (L/N_PROC_X)
#define LY (L)
#define LZ (L)
#define VOLUME (T*LX*LY*LZ)

#define DUM_DERI 6
#define DUM_SOLVER (DUM_DERI+4)
#define DUM_MATRIX (DUM_SOLVER+6) 
/* if you want to include bicgstabell */
/* #define DUM_MATRIX DUM_SOLVER+11 */
#define NO_OF_SPINORFIELDS (DUM_MATRIX+2)
/* For benchmark set the following: */
/* #define NO_OF_SPINORFIELDS 100  */

#if (defined PARALLELT && !defined PARALLELXT)
#define RAND (2*LX*LY*LZ)
#define VOLUMEPLUSRAND ((T+2)*LX*LY*LZ)
#elif defined PARALLELXT
#define RAND (2*LY*LZ*(LX+T))
/* Note that VOLUMEPLUSRAND not equal to VOLUME+RAND in this case */
#define VOLUMEPLUSRAND (LY*LZ*(T+2)*(LX+2))
#else
#define RAND 0
#define VOLUMEPLUSRAND (VOLUME)
#endif

/* Here you can define antiperiodic  */
/* boundary conditions with e.g.     */
/* #define X1 1.  (in time)          */
#define X1 0.
#define X2 0.
#define X3 0.
#define EPS_SQ0  1.0e-20
#define EPS_SQ1  1.0e-7
#define EPS_SQ2  1.0e-5
#define EPS_SQ3  1.0e-3
#define tiny_t  1.0e-20

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

#if ((defined SSE)||(defined SSE2)||(defined SSE3))
#include "sse.h"
#else
  #define ALIGN
#endif

/* translates from lexicagraphic order to even/odd order */
EXTERN int * g_lexic2eo;
/* translates from even/odd orderto lexicograhic order  */    
EXTERN int * g_eo2lexic;
EXTERN int * g_lexic2eosub;

EXTERN int **** g_ipt;
EXTERN int ** g_iup;
EXTERN int ** g_idn;

EXTERN spinor ** spinor_field;
/* EXTERN spinor ** spinor_field; */
EXTERN su3 ** g_gauge_field;
EXTERN su3 ** g_gauge_field_back;
/* This is dirty, but dow not allocate memory */
/* if no clover is used. */
EXTERN su3adj ** moment;
EXTERN su3adj ** df0;
EXTERN su3adj ** ddummy;
#ifdef CLOVER
EXTERN su3adj dclover[VOLUMEPLUSRAND][4] ALIGN;
EXTERN su3 sw[VOLUME][3][2] ALIGN;
EXTERN su3 sw_inv[VOLUME][3][2] ALIGN;
EXTERN su3 swp[VOLUME][4] ALIGN;
EXTERN su3 swm[VOLUME][4] ALIGN;
#else
EXTERN su3adj dclover[1][1] ALIGN;
EXTERN su3 sw[1][1][1] ALIGN;
EXTERN su3 sw_inv[1][1][1] ALIGN;
EXTERN su3 swp[1][1] ALIGN;
EXTERN su3 swm[1][1] ALIGN;
#endif
EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double g_kappa, g_c_sw, g_ka_csw_8, g_beta;
EXTERN double g_rgi_C0, g_rgi_C1;
EXTERN double g_mu, g_mu1, g_mu2, g_mu3;
EXTERN int g_use_clover_flag, g_nr_of_psf;

/* MPI information */
EXTERN int g_proc_id, g_nproc, g_stdio_proc, g_nproc_t, g_nproc_x, g_cart_id;
EXTERN int g_proc_coords[3];
#ifdef MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
EXTERN MPI_Comm g_cart_grid;

/* the next neighbours for MPI */
EXTERN int g_nb_x_up, g_nb_x_dn;
EXTERN int g_nb_t_up, g_nb_t_dn;


#endif

#undef EXTERN
/* #undef ALIGN */

#ifdef MAIN_PROGRAM
static char const g_rcsid[] = "$Id$";
#endif

#endif

