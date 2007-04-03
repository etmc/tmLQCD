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

void jderrorhandler(const int i, char * message);

extern void jdher(int n, int lda, double tau, double tol, 
		  int kmax, int jmax, int jmin, int itmax,
		  int blksize, int blkwise, 
		  int V0dim, complex *V0, 
		  int solver_flag, 
		  int linitmax, double eps_tr, double toldecay,
		  int verbosity,
		  int *k_conv, complex *Q, double *lambda, int *it,
		  int maxmin, int shift_mode,
		  matrix_mult A_psi);

#endif
