#ifndef _GRAM_SCHMIDT_BI_H
#define _GRAM_SCHMIDT_BI_H
#include "complex.h"

void IteratedClassicalGS_bi_old(complex v[], double *vnrm, int n, int m, complex A[], complex work1[]);
void IteratedClassicalGS_bi(complex v[], double *vnrm, int n, int m, complex A[], 
			 complex work1[], int lda) ;

void ModifiedGS_bi_old(complex v[], int n, int m, complex A[]);
void ModifiedGS_bi(complex v[], int n, int m, complex A[], int lda);

#endif
