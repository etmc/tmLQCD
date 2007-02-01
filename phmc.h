/* $Id$ */

#ifndef _PHMC_H
#define _PHMC_H

#ifndef MAIN_PROGRAM
extern double phmc_Cpol;
extern double phmc_cheb_evmin, phmc_cheb_evmax;
extern double phmc_invmaxev;
extern complex * phmc_roo;
extern int phmc_dop_n_cheby;
extern double * phmc_dop_cheby_coef;
extern int phmc_ptilde_n_cheby;
extern double * phmc_ptilde_cheby_coef;
extern double phmc_stilde_low, phmc_stilde_max;
#else
double phmc_Cpol;
double phmc_cheb_evmin, phmc_cheb_evmax;
double phmc_invmaxev;
complex * phmc_roo;
int phmc_dop_n_cheby;
double * phmc_dop_cheby_coef;
int phmc_ptilde_n_cheby;
double * phmc_ptilde_cheby_coef;
double phmc_stilde_low, phmc_stilde_max;
#endif


#endif
