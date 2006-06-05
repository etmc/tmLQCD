/* $Id$ */
#ifndef _PTILDE_ND_H
#define _PTILDE_ND_H


double func_tilde(double u, double exponent);

void Ptilde_cheb_coefs(double a, double b, double dd[], int n, double exponent);

void Poly_tilde_ND(spinor *R_s, spinor *R_c, double *dd, int n, spinor *S_s, spinor *S_c);

void degree_of_Ptilde();

#endif
