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

#define T  2
#define L  4

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
#define DUM_SOLVER DUM_DERI+4
#define DUM_MATRIX DUM_SOLVER+5 
/* if you want to include bicgstabell */
/* #define DUM_MATRIX DUM_SOLVER+11 */
#define NO_OF_SPINORFIELDS DUM_MATRIX+2

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

#if ((defined SSE)||(defined SSE2))
  #if defined P4
    #define ALIGN_BASE 0x3f
    #define ALIGN __attribute__ ((aligned (64)))
  #else
    #define ALIGN_BASE 0x1f
    #define ALIGN __attribute__ ((aligned (32)))
  #endif
#else
  #define ALIGN
#endif

/* translates from lexicagraphic order to even/odd order */
EXTERN int g_lexic2eo[VOLUMEPLUSRAND] ALIGN;
/* translates from even/odd orderto lexicograhic order  */    
EXTERN int g_eo2lexic[VOLUMEPLUSRAND] ALIGN;
EXTERN int g_lexic2eosub[VOLUMEPLUSRAND] ALIGN;
EXTERN int g_ipt[T+2][LX+2][LY][LZ] ALIGN;
EXTERN int g_iup[VOLUMEPLUSRAND][4] ALIGN;   
EXTERN int g_idn[VOLUMEPLUSRAND][4] ALIGN;   
EXTERN spinor spinor_field[NO_OF_SPINORFIELDS][(VOLUMEPLUSRAND)/2] ALIGN;
/* EXTERN spinor ** spinor_field; */
EXTERN su3 g_gauge_field[VOLUMEPLUSRAND][4] ALIGN;
EXTERN su3adj moment[VOLUME][4] ALIGN;
EXTERN su3adj df0[VOLUMEPLUSRAND][4] ALIGN;
EXTERN su3adj dclover[VOLUMEPLUSRAND][4] ALIGN;
EXTERN su3adj ddummy[VOLUMEPLUSRAND][4] ALIGN;
EXTERN su3 sw[VOLUME][3][2] ALIGN;
EXTERN su3 sw_inv[VOLUME][3][2] ALIGN;
EXTERN su3 swp[VOLUME][4] ALIGN;
EXTERN su3 swm[VOLUME][4] ALIGN;
EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double g_kappa, g_c_sw, g_ka_csw_8, g_beta;
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

#endif
