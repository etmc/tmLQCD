/* $Id$ */
#ifndef _CHEBYSHEV_POLYNOMIAL_H
#define _CHEBYSHEV_POLYNOMIAL_H

double func(double u);
void chebyshev_polynomial(double a, double b, double c[], int n);

void QdaggerQ_onequarter(spinor *R_s, spinor *R_c, double *c, int n, spinor *S_s, spinor *S_c, double minev);

void degree_of_polynomial();

#endif
