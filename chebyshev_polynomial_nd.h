/* $Id$ */
#ifndef _CHEBYSHEV_POLYNOMIAL_ND_H
#define _CHEBYSHEV_POLYNOMIAL_ND_H


double func(double u, double exponent);

void chebyshev_coefs(double a, double b, double c[], int n, double exponent);

void QdaggerQ_poly(spinor *R_s, spinor *R_c, double *c, int n, spinor *S_s, spinor *S_c);

void degree_of_polynomial_nd();

#endif
