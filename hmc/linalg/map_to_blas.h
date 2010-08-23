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
/* $Id$ */

/*
 * Mapping from our linalg routines
 * to blas routines
 *
 * Carsten Urbach May 2003
 * urbach@ifh.de
 */

#ifndef MAP_TO_BLAS_H
#define MAP_TO_BLAS_H

#ifdef _USE_BLAS

#ifdef XLC
/*#include <essl.h>*/
#include "su3/complex.h"
complex zdotc(int, complex*, int, complex*, int);
void zaxpy(int, complex ,complex* ,int ,complex* ,int);
void zcopy(int, complex*, int, complex*, int);
#define assign_add_mul(A,B,C) zaxpy(12*VOLUME,C,(complex*)B,1,(complex*)A,1)
#define scalar_prod(A,B) zdotc(12*VOLUME,(complex*)A,1,(complex*)B,1)
#define assign(A,B) zcopy(12*VOLUME,(complex*)B,1,(complex*)A,1)

#else

#define assign_add_mul(A,B,C) zaxpy(12*VOLUME,C,B,1,A,1)
#define scalar_prod(A,B) zdotc(12*VOLUME,A,1,B,1)
#define assign(A,B) zcopy(12*VOLUME,B,1,A,1)
#endif

#endif
#endif
