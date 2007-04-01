/* $Id$ */
#ifndef _EIGENVALUES_H
#define _EIGENVALUES_H

extern spinor * eigenvectors;
extern double * eigenvls;
extern int eigenvalues_for_cg_computed;
void eigenvalues(int * nev, const int max_iterations, const double prec);

#endif
