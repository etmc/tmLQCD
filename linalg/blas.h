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

#ifndef _BLAS_H
#define _BLAS_H

#include <complex.h>
#include "linalg/fortran.h"

#if defined CRAY || defined HITACHI
/* On the CRAY is all different, of course... */
#include"fortran.h"
#define zgemm ZGEMM
#define zgemv ZGEMV
#define ddot DDOT
#define zdotc ZDOTC
#define daxpy DAXPY
#define dnrm2 DNRM2
#define znrm2 ZNRM2
#define zaxpy ZAXPY
#define dcopy DCOPY
#define dscal DSCAL
#define dgemv DGEMV
#define dgemm DGEMM
extern double _FT(dasum);
extern double _FT(ddot)();
extern void _FT(zdotc)();
extern double _FT(dnrm2)();
extern double _FT(znrm2)();
extern int _FT(idamax)();
extern void _FT(daxpy)();
extern void _FT(zaxpy)();
extern void _FT(dcopy)();
extern void _FT(dscal)();
extern void _FT(dgemv)();
extern void _FT(zgemv)();
extern void _FT(dgemm)();
extern void _FT(zgemm)();
#else

/* BLAS-1 functions */
extern double _FT(dasum)(int* n, double x[], int* incx);
extern double _FT(ddot)(int* n, double x[], int* incx, double y[],
        int* incy);
extern void _FT(zdotc)(int* n, _Complex double x[], int* incx, _Complex double y[],
        int* incy);
extern double _FT(dnrm2)(int* n, double x[], int* incx);
extern double _FT(znrm2)(int* n, _Complex double x[], int* incx);
extern int _FT(idamax)(int* n, double x[], int* incx);

/* BLAS-1 subroutines */
extern void _FT(daxpy)(int* n, double* a, double x[], int* incx,
        double y[], int* incy);
extern void _FT(zaxpy)(int* n, _Complex double* a, _Complex double x[], int* incx,
        _Complex double y[], int* incy);
extern void _FT(dcopy)(int* n, double x[], int* incx, double y[],
        int* incy);
extern void _FT(dscal)(int* n, double* a, double x[], int* incx);

/* BLAS-2 subroutines */
extern void _FT(dgemv)(char* trans, int* m, int* n, double* alpha,
        double a[], int* lda, double x[], int* incx, double* beta,
        double y[], int* incy, int len_trans);
extern void _FT(zgemv)(char* trans, int* m, int* n, _Complex double* alpha,
        _Complex double a[], int* lda, _Complex double x[], int* incx, _Complex double* beta,
        _Complex double y[], int* incy, int len_trans);

/* BLAS-3 subroutines */
extern void _FT(dgemm)(char* transa, char* transb, int* m, int* n, int* k,
        double* alpha, double a[], int* lda, double b[], int* ldb,
        double* beta, double c[], int* ldc, int len_transa, int len_transb);
extern void _FT(zgemm)(char* transa, char* transb, int* m, int* n, int* k,
        _Complex double* alpha, _Complex double a[], int* lda, _Complex double b[], int* ldb,
        _Complex double* beta, _Complex double c[], int* ldc, int len_transa, int len_transb);
#endif

#endif
