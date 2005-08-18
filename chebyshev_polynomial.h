/* $Id$ */
#ifndef _CHEBYSHEV_POLYNOMIAL_H
#define _CHEBYSHEV_POLYNOMIAL_H

extern double cheb_evmin, cheb_evmax;
extern int dop_n_cheby;
extern double * dop_cheby_coef;


double func(double u, double exponent);
void chebyshev_polynomial(double a, double b, double c[], int n, double exponent);

void QdaggerQ_power(spinor *R_s, spinor *R_c, double *c, int n, spinor *S_s, spinor *S_c);

void degree_of_polynomial();

#endif
