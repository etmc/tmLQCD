/* $Id$ */
#ifndef _EIGENVALUES_H
#define _EIGENVALUES_H

extern spinor * eigenvectors;
extern double * eigenvls;
extern int eigenvalues_for_cg_computed;
double eigenvalues(int * nr_of_eigenvalues, const int max_iterations, 
		   const double precision, const int maxmin,
		   const int readwrite, const int nstore);

#endif
