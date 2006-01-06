/* $Id$ */
/*******************************************************************************
*
* File linalg.c
*
* Linear combination and norm of spinor fields
*
* The externally accessible function are
*
*   void assign_add_mul_r(double c,int k,int l)
*     Adds c*psi[k] to psi[l]
*
*   double square_norm(int k)
*     Returns the square norm of psi[k]
*
* Based on
* Version: 2.0
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 21.03.2001
* Present version: M.Hasenbusch Thu Oct 25 10:58:45 METDST 2001
* Extended and adapted for even-odd precondtioning
*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "linalg_eo.h"

/* j output k input , l input */




/* j, k input, l output */












/* j output k input , l input */



