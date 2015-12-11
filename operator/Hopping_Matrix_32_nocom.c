/***********************************************************************
 * Copyright (C) 2013 Florian Burger
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

// work-around for missing single precision implementation of inline SSE
#ifdef SSE
#define REDEFSSE
#undef SSE
#endif

#ifdef SSE2
#define REDEFSSE2
#undef SSE2
#endif

#ifdef SSE3
#define REDEFSSE3
#undef SSE3
#endif

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "xchange/xchange.h"
#include "su3.h"
#include "sse.h"
#include "boundary.h"
#include "operator/Hopping_Matrix_32.h"

#define Hopping_Matrix_32 Hopping_Matrix_32_nocom
#define _NO_COMM 1

#include "Hopping_Matrix_32.c"
