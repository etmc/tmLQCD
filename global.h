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

#define T  4
#define L  8

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

#ifdef PARALLELT
#define RAND (2*LX*LY*LZ)
#define VOLUMEPLUSRAND ((T+2)*LX*LY*LZ)
#elif defined PARALLELXT
#define RAND (2*LY*LZ*(LX+T))
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
#if ((defined SSE)||(defined SSE2))
  #if defined P4
    #define ALIGN __attribute__ ((aligned (64)))
  #else 
    #define ALIGN __attribute__ ((aligned (32)))
  #endif
#else
  #define ALIGN
#endif
#else
  #define EXTERN extern
  #define ALIGN
#endif

/* translates from ix to ieven/odd  */
EXTERN int trans1[VOLUME+RAND] ALIGN;
/* translates from ieven/iodd to ix */    
EXTERN int trans2[VOLUME+RAND] ALIGN;   
EXTERN int xeven[VOLUME+RAND] ALIGN;   
EXTERN int g_ipt[T+2][LX+2][LY][LZ] ALIGN;
EXTERN int g_iup[VOLUME+RAND][4] ALIGN;   
EXTERN int g_idn[VOLUME+RAND][4] ALIGN;   
EXTERN spinor spinor_field[NO_OF_SPINORFIELDS][(VOLUME+RAND)/2] ALIGN;
EXTERN su3 g_gauge_field[VOLUME+RAND][4] ALIGN;
EXTERN su3adj moment[VOLUME][4] ALIGN;
EXTERN su3adj df0[VOLUME+RAND][4] ALIGN;
EXTERN su3adj dclover[VOLUME+RAND][4] ALIGN;
EXTERN su3adj ddummy[VOLUME+RAND][4] ALIGN;
EXTERN su3 sw[VOLUME][3][2] ALIGN;
EXTERN su3 sw_inv[VOLUME][3][2] ALIGN;
EXTERN su3 swp[VOLUME][4] ALIGN;
EXTERN su3 swm[VOLUME][4] ALIGN;
EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double g_kappa, g_c_sw, g_ka_csw_8, g_beta;
EXTERN double g_mu, g_mu1, g_mu2, g_mu3;
EXTERN int g_use_clover_flag, g_nr_of_psf;

/* MPI information */
EXTERN int g_proc_id, g_nproc, g_stdio_proc, g_nproc_t;
EXTERN int g_proc_coords[3];
#ifdef MPI
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;
EXTERN MPI_Comm g_cart_grid;

/* the next neighbours for MPI */
EXTERN int g_nb_x_up, g_nb_x_dn;
EXTERN int g_nb_t_up, g_nb_t_dn;

/* Datatypes for the data exchange */
EXTERN MPI_Datatype g_gauge_point, g_gauge_time_slice_cont;
EXTERN MPI_Datatype g_gauge_time_slice_split;
EXTERN MPI_Datatype g_deri_point, g_deri_time_slice_cont;
EXTERN MPI_Datatype g_deri_time_slice_split;
EXTERN MPI_Datatype g_field_point, g_field_time_slice_cont;
EXTERN MPI_Datatype g_gauge_x_slice_cont;
EXTERN MPI_Datatype g_gauge_yz_subslice;
EXTERN MPI_Datatype g_gauge_x_slice_gath;
EXTERN MPI_Datatype g_field_x_slice_cont;
EXTERN MPI_Datatype g_field_yz_subslice;
EXTERN MPI_Datatype g_field_x_slice_gath;
EXTERN MPI_Datatype g_deri_x_slice_cont;
EXTERN MPI_Datatype g_deri_yz_subslice;
EXTERN MPI_Datatype g_deri_x_slice_gath;

#endif

#undef EXTERN
/* #undef ALIGN */

#endif
