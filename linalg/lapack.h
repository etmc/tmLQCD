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

#ifndef _LAPACK_H
#define _LAPACK_H

#include <complex.h>
#include "linalg/fortran.h"

#if defined CRAY || defined HITACHI
#define zgels CGELS
#define zgesv CGESV
#define zgeevx CGEEVX
#define dsyev CSYEV
#define zheev CHEEV
#define dgetrs DGETRS
#define dgetrf DGETRF
#define dlarnv DLARNV
#define zlarnv CLARNV
#define dsyevx DSYEVX
#define zlacpy CLACPY
#define dlacpy DLACPY
#define dlaset DLASET
#define zlaset CLASET 
#define dlamch DLAMCH
#define ilaenv ILAENV
#define zlapcy CLAPCY
#define zgetrf CGETRF
#define zgetrs CGETRS
#define zgeqrf ZGEQRF
#define zunmqr ZUNMQR

extern void _FT(zgels)();
extern void _FT(zgesv)();
extern void _FT(zgeevx)();
extern void _FT(dsyev)();
extern void _FT(zheev)();
extern void _FT(dgetrs)();
extern void _FT(dgetrf)();
extern void _FT(dlarnv)();
extern void _FT(zlarnv)();
extern void _FT(dsyevx)();
extern void _FT(zlacpy)();
extern void _FT(dlaset)();
extern double _FT(dlamch)();
extern int _FT(ilaenv)();
extern void _FT(zgetrf)();
extern void _FT(zgetrs)();
extern void _FT(zgeqrf)();
extern void _FT(zunmqr)();

#else

void _FT(zgels)(char* transa, int* M, int* N, int* NRHS, _Complex double a[], int* lda, 
		 _Complex double b[], int* ldb, _Complex double work[], int* lwork, int* info, int len_transa);

void _FT(zgesv)(int* n, int* nrhs, _Complex double a[], int* lda,
		 int ipivot[], _Complex double b[], int* ldb, int *info);

extern void _FT(zgeevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
			 int* N, _Complex double A[], int* lda, _Complex double W[], _Complex double vl[], 
			 int* ldvl, _Complex double vr[], int* ldvr, int* ilo, int* ihi,
			 double scale[], double* abnrm, double rcone[], double rconv[],
			 _Complex double work[], int* lwork, double work2[], int* info, 
			 int len_balanc, int len_jobvl, int len_jobvr, int len_sense);

extern void _FT(dsyev)(char* jobz, char* uplo, int* n, double a[],
        int* lda, double w[], double work[], int* lwork, int* info,
        int len_jobz, int len_uplo);
extern void _FT(zheev)(char* jobz, char* uplo, int* n, _Complex double a[],
        int* lda, double w[], _Complex double work[], int* lwork, double* rwork, int* info,int len_jobz, int len_uplo);

extern void _FT(dgetrs)(char* trans, int* n, int* nrhs, double a[],
        int* lda, int ipiv[], double b[], int* ldb, int* info,
        int len_trans);
extern void _FT(dgetrf)(int* m, int* n, double a[], int* lda, int ipiv[],
        int* info);

extern void _FT(zgetrs)(char* trans, int* n, int* nrhs, _Complex double a[],
        int* lda, int ipiv[], _Complex double b[], int* ldb, int* info,
        int len_trans);
extern void _FT(zgetrf)(int* m, int* n, _Complex double a[], int* lda, int ipiv[],
        int* info);

extern void _FT(zhetrs)(char* uplo, int* n, int* nrhs, _Complex double a[],
			 int* lda, int ipiv[], _Complex double b[], int* ldb, int* info,
			 int len_uplo);
extern void _FT(zhetrf)(char* uplo, int* n, _Complex double a[], int* lda, int ipiv[],
			 _Complex double work[], int * lwork, int* info, int len_uplo);

extern void _FT(dlarnv)(int *IDIST, int *ISEED, int *N, double *X);
extern void _FT(zlarnv)(int *IDIST, int *ISEED, int *N, _Complex double *X);

extern void _FT(dsyevx)(char* jobz, char* range, char* uplo, int* n,
        double a[], int* lda, double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], double z[], int* ldz,
        double work[], int* lwork, int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void _FT(zlacpy)(char *UPLO, int *M, int *N, _Complex double *A, int *LDA, 
			_Complex double *B, int *LDB, int len_uplo);
extern void _FT(dlacpy)(char *UPLO, int *M, int *N, double *A, int *LDA, 
			double *B, int *LDB, int len_uplo);

extern void _FT(dlaset)(char *UPLO, int *M, int *N, double *ALPHA, 
			double *BETA, double *A, int *LDA, int len_uplo );
extern void _FT(zlaset)(char *UPLO, int *M, int *N, _Complex double *ALPHA, 
			_Complex double *BETA, _Complex double *A, 
                        int *LDA, int len_uplo );

extern double _FT(dlamch)(char* name, int len_name);

extern int _FT(ilaenv)(int *ISPEC, char *NAME, char *OPTS, int *N1, 
		       int *N2, int *N3, int *N4, int len_name, int len_opts);

extern void _FT(zgeqrf)(int *M, int *N, _Complex double *A, int *LDA, _Complex double *TAU,
                         _Complex double  *WORK, int *LWORK, int *INFO);


extern void _FT(zunmqr)(char *SIDE, char *TRANS, int *M, int *N, int *K,
                         _Complex double  *A, int *LDA, _Complex double  *TAU, _Complex double  *C,
                         int *LDC, _Complex double  *WORK, int *LWORK, int *INFO);
#endif

#endif
