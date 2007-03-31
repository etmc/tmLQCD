/* $Id$ */
#ifndef _EIGENVALUES_BI_H
#define _EIGENVALUES_BI_H

extern bispinor * eigenvectors;
extern double * eigenvls;

double eigenvalues_bi(int * nev, const int operator_flag, const int max_iterations, 
		      const double prec, const int maxmin);

#endif
