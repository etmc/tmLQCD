/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
		  int V0dim, _Complex double *V0, 
		  int solver_flag, 
		  int linitmax, double eps_tr, double toldecay,
		  int verbosity,
		  int *k_conv, _Complex double *Q, double *lambda, int *it,
		  int maxmin, int shift_mode,
		  matrix_mult A_psi);

#endif
