/***********************************************************************
 *
 * Copyright (C) 2008 Albert Deuzeman, Siebren Reker, Carsten Urbach
 *               2010 Claude Tadonki, Carsten Urbach
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


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <complex.h>
#include "block.h"
#include "little_field_gather.h"

enum{
  NONE = 0,
  T_UP = 1,
  T_DN = 2,
  X_UP = 3,
  X_DN = 4,
  Y_UP = 5,
  Y_DN = 6,
  Z_UP = 7,
  Z_DN = 8
} Direction;


#ifdef MPI
MPI_Request lrequests[16];
MPI_Status lstatus[16];
int waitcount = 0;
#endif


#define _PSWITCH(s) s 
#define _C_TYPE _Complex double
#define _MPI_C_TYPE MPI_DOUBLE_COMPLEX

#include"little_field_gather_body.c"

#undef _PSWITCH
#undef _C_TYPE
#undef _MPI_C_TYPE

#define _PSWITCH(s) s ## 32
#define _C_TYPE _Complex float
#define _MPI_C_TYPE MPI_COMPLEX

#include"little_field_gather_body.c"

#undef _PSWITCH
#undef _C_TYPE
#undef _MPI_C_TYPE
