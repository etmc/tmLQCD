/* $Id$ */

#ifndef _JDHER_H
#define _JDHER_H

#ifndef JD_MAXIMAL
#define JD_MAXIMAL 1
#endif
#ifndef JD_MINIMAL
#define JD_MINIMAL 0
#endif

#include <stdlib.h>
#include <stdio.h>
#include "su3.h"
#include "solver/solver.h"

#define errorhandler(a, b) \
        {printf("Error %d %s!\n", a, b); \
        fflush(stdout); \
        exit(6);}

extern void jdher(int n, double tau, double jdtol, 
		  int kmax, int jmax, int jmin, int itmax,
		  int blksize, int blkwise, 
		  int V0dim, complex *V0, 
		  int linsolver,  
		  int linitmax, double eps_tr, double toldecay,
		  int clvl,
		  int *k_conv, complex *Q, double *lambda, int *it,
		  int maxmin, const int shift_mode,
		  matrix_mult domatveca);

#endif
