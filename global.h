
/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
* Version: 1.1
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 16.03.2001
*
* Adapted for the HMC-Program by M. Hasenbusch 2002
*
*******************************************************************************/
#include "mpi.h"

#define T  4
#define L  8
#define DUM_DERI 6
#define DUM_SOLVER DUM_DERI+4
#define DUM_MATRIX DUM_SOLVER+5
#define NO_OF_SPINORFIELDS DUM_MATRIX+2
#define VOLUME (T*L*L*L)
#define RAND (2*L*L*L)
#define SLICE (L*L*L/2)
/* #define beta  5.1 */
/* #define kappa 0.12 */
/* #define ka_csw_8 0.03058 */
/* #define kappa 0.1 */
/* #define ka_csw_8 0.01875   */
/* #define ka_csw_8 0.0  */
/* #define X0 1. */
#define X0 0.
#define X1 0.
#define X2 0.
#define X3 0.
#define ITER_MAX 5000
/* #define ITER_MAX 1 */
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

EXTERN int trans1[VOLUME+RAND] ALIGN;   
EXTERN int trans2[VOLUME+RAND] ALIGN;   
EXTERN int xeven[VOLUME+RAND] ALIGN;   
EXTERN int ipt[T+2][L][L][L] ALIGN;
EXTERN int iup[VOLUME+RAND][4] ALIGN;   
EXTERN int idn[VOLUME+RAND][4] ALIGN;   
EXTERN spinor spinor_field[NO_OF_SPINORFIELDS][(VOLUME+RAND)/2] ALIGN;
EXTERN su3 gauge_field[VOLUME+RAND][4] ALIGN;
EXTERN su3adj moment[VOLUME][4] ALIGN;
EXTERN su3adj df0[VOLUME+RAND][4] ALIGN;
EXTERN su3adj dclover[VOLUME+RAND][4] ALIGN;
EXTERN su3adj ddummy[VOLUME+RAND][4] ALIGN;
EXTERN su3 sw[VOLUME][3][2] ALIGN;
EXTERN su3 sw_inv[VOLUME][3][2] ALIGN;
EXTERN su3 swp[VOLUME][4] ALIGN;
EXTERN su3 swm[VOLUME][4] ALIGN;
EXTERN int count00,count01,count10,count11,count20,count21;
EXTERN double kappa, c_sw, ka_csw_8, beta;

/* MPI information */
EXTERN int myid, numprocs;
EXTERN MPI_Status status;
EXTERN MPI_Request req1,req2,req3,req4;

/* file for solver-information */
FILE *fp7;

#undef EXTERN
#undef ALIGN
