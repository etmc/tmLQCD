/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2010 claude Tadonki
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
#include "tmlqcd_config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef TM_USE_OMP
#include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
// needed for the dummy function little_mg_precon ...
#include "block.h"
#include "dfl_projector.h"
#include "gcr4complex.h"


#define _PSWITCH(s) s 
#define _PTSWITCH(s) s 
#define _C_TYPE _Complex double
#define _F_TYPE double

#include "gcr4complex_body.c"

#undef _PSWITCH
#undef _PTSWITCH
#undef _C_TYPE
#undef _F_TYPE

#define _PSWITCH(s) s ## _32
#define _PTSWITCH(s) s ## 32
#define _C_TYPE _Complex float
#define _F_TYPE float

#include "gcr4complex_body.c"

#undef _PSWITCH
#undef _PTSWITCH
#undef _C_TYPE
#undef _F_TYPE
