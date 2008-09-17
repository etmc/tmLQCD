/* $Id$*/

#ifndef LU_SOLVE_H
#define LU_SOLVE_H

/* Solve M a = b by LU decomposition with partial pivoting */
void LUSolve( const int Nvec, complex * M, const int ldM, complex * b);

void LUInvert( const int Nvec, complex * const M, const int ldM);

#endif
