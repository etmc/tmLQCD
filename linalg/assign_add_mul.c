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
/*******************************************************************************
 *
 * File assign_add_mul.c 
 *
 *   void assign_add_mul(spinor * const P, spinor * const Q, const complex c)
 *     (*P) = (*P) + c(*Q)        c is a complex constant
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "assign_add_mul.h"


#define _C_TYPE _Complex double
#define _PSWITCH(s) s
#define _PTSWITCH(s) s

#include "assign_add_mul_body.c"

#undef _C_TYPE
#undef _PSWITCH
#undef _PTSWITCH

#define _C_TYPE _Complex float
#define _PSWITCH(s) s ## _32
#define _PTSWITCH(s) s ## 32

#include "assign_add_mul_body.c"

#undef _C_TYPE
#undef _PSWITCH
#undef _PTSWITCH

