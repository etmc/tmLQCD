/* $Id$ */

#ifndef _LAPACK_H
#define _LAPACK_H

#include "linalg/fortran.h"

#if defined CRAY || defined HITACHI
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
#else

extern void _FT(dsyev) (char* jobz, char* uplo, int* n, double a[],
        int* lda, double w[], double work[], int* lwork, int* info,
        int len_jobz, int len_uplo);
extern void _FT(zheev) (char* jobz, char* uplo, int* n, complex a[],
        int* lda, double w[], complex work[], int* lwork, double* rwork, int* info,
        int len_jobz, int len_uplo);

extern void _FT(dgetrs) (char* trans, int* n, int* nrhs, double a[],
        int* lda, int ipiv[], double b[], int* ldb, int* info,
        int len_trans);
extern void _FT(dgetrf) (int* m, int* n, double a[], int* lda, int ipiv[],
        int* info);

extern void _FT(zgetrs) (char* trans, int* n, int* nrhs, complex a[],
        int* lda, int ipiv[], complex b[], int* ldb, int* info,
        int len_trans);
extern void _FT(zgetrf) (int* m, int* n, complex a[], int* lda, int ipiv[],
        int* info);

extern void _FT(zhetrs) (char* uplo, int* n, int* nrhs, complex a[],
			 int* lda, int ipiv[], complex b[], int* ldb, int* info,
			 int len_uplo);
extern void _FT(zhetrf) (char* uplo, int* n, complex a[], int* lda, int ipiv[],
			 complex work[], int * lwork, int* info, int len_uplo);

extern void _FT(dlarnv) (int *IDIST, int *ISEED, int *N, double *X);
extern void _FT(zlarnv) (int *IDIST, int *ISEED, int *N, complex *X);

extern void _FT(dsyevx) (char* jobz, char* range, char* uplo, int* n,
        double a[], int* lda, double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], double z[], int* ldz,
        double work[], int* lwork, int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void _FT(zlacpy)(char *UPLO, int *M, int *N, complex *A, int *LDA, 
			complex *B, int *LDB, int len_uplo);
extern void _FT(dlacpy)(char *UPLO, int *M, int *N, double *A, int *LDA, 
			double *B, int *LDB, int len_uplo);

extern void _FT(dlaset)(char *UPLO, int *M, int *N, double *ALPHA, 
			double *BETA, double *A, int *LDA, int len_uplo );
extern void _FT(zlaset)(char *UPLO, int *M, int *N, complex *ALPHA, 
			complex *BETA, complex *A, int *LDA, int len_uplo );

extern double _FT(dlamch) (char* name, int len_name);

extern int _FT(ilaenv)(int *ISPEC, char *NAME, char *OPTS, int *N1, 
		       int *N2, int *N3, int *N4, int len_name, int len_opts);

#endif

#endif
