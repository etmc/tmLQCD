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
