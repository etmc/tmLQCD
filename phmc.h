/* $Id$ */

#ifndef _PHMC_H
#define _PHMC_H

#ifndef MAIN_PROGRAM
/* the normalisation constant appearing in the product representation of */
/* the polynomial */
extern double phmc_Cpol;
/* maximal and minimal eigenvalue of the ND operator */
extern double phmc_cheb_evmin, phmc_cheb_evmax;
/* inverse maximal EV, needed for normalisation */
extern double phmc_invmaxev;
/* These are the roots */
extern complex * phmc_roo;
/* degree and coefs of P */
extern int phmc_dop_n_cheby;
extern double * phmc_dop_cheby_coef;
/* degree of coefs \tilde P */
extern int phmc_ptilde_n_cheby;
extern double * phmc_ptilde_cheby_coef;
#else
double phmc_Cpol;
double phmc_cheb_evmin, phmc_cheb_evmax;
double phmc_invmaxev;
complex * phmc_roo;
int phmc_dop_n_cheby;
double * phmc_dop_cheby_coef;
int phmc_ptilde_n_cheby;
double * phmc_ptilde_cheby_coef;
#endif


#endif
