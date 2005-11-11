/* $Id$ */
#ifndef _MAX_EIGENVALUES_BI_H
#define _MAX_EIGENVALUES_BI_H

extern bispinor * max_evs;
extern double * max_evls;
/*
extern int eigenvalues_for_cg_computed;
*/
void max_eigenvalues(int * nev, const int operator_flag, const int max_iterations, const double prec);

#endif
