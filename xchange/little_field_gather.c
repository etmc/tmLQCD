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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef TM_USE_MPI
# include <mpi.h>
#endif
#include "global.h"
#include <complex.h>
#include "block.h"
#include "little_field_gather.h"

enum{
  T_UP = 0,
  T_DN = 1,
  X_UP = 2,
  X_DN = 3,
  Y_UP = 4,
  Y_DN = 5,
  Z_UP = 6,
  Z_DN = 7
} GatherDirection;


#ifdef TM_USE_MPI
MPI_Request lrequests[16];
MPI_Status lstatus[16];
int waitcount = 0;
#endif


#define _PSWITCH(s) s 
#define _PTSWITCH(s) s 
#define _C_TYPE _Complex double
#define _MPI_C_TYPE MPI_DOUBLE_COMPLEX

#include"little_field_gather_body.c"

#undef _PSWITCH
#undef _PTSWITCH
#undef _C_TYPE
#undef _MPI_C_TYPE

#define _PSWITCH(s) s ## _32
#define _PTSWITCH(s) s ## 32
#define _C_TYPE _Complex float
#define _MPI_C_TYPE MPI_COMPLEX

#include"little_field_gather_body.c"

#undef _PSWITCH
#undef _PTSWITCH
#undef _C_TYPE
#undef _MPI_C_TYPE
