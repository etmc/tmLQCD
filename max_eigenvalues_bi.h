/* $Id$ */
#ifndef _MAX_EIGENVALUES_BI_H
#define _MAX_EIGENVALUES_BI_H

extern bispinor * max_evs;
extern double * max_evls;
/*
extern int eigenvalues_for_cg_computed;
*/
double max_eigenvalues_bi(int * nev, const int operator_flag, const int max_iterations, const double prec);

#endif
