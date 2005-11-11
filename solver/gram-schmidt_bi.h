#ifndef _GRAM_SCHMIDT_BI_H
#define _GRAM_SCHMIDT_BI_H
#include "complex.h"

void IteratedClassicalGS(complex v[], double *vnrm, int n, int m, complex A[], complex work1[]);
void pIteratedClassicalGS(complex v[], double *vnrm, int n, int m, complex A[], 
			 complex work1[], int lda) ;

void ModifiedGS(complex v[], int n, int m, complex A[]);
void pModifiedGS(complex v[], int n, int m, complex A[], int lda);

#endif
