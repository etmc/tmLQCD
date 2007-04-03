#ifndef _GRAM_SCHMIDT_H
#define _GRAM_SCHMIDT_H
#include "complex.h"

void IteratedClassicalGS_old(complex v[], double *vnrm, int n, int m, complex A[], complex work1[]);
void IteratedClassicalGS(complex v[], double *vnrm, int n, int m, complex A[], 
			 complex work1[], int lda) ;

void ModifiedGS_old(complex v[], int n, int m, complex A[]);
void ModifiedGS(complex v[], int n, int m, complex A[], int lda);

#endif
