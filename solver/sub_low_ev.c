/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
/******************************************************************** 
 *
 * void sub_low_ev(spinor *S, spinor *P)
 * makes
 * |S> = |P> - Sum_{1}^{nev-1} <eigen_i|P>*|eigen_i> 
 *  
 * where |eigen_i> is the i-th lowest eigenvectors of Q 
 * 
 *
 * void addproj_q_invsqrt(spinor *Q, spinor *P)
 * makes
 * |Q'> = |Q> + Sum_{1}^{nev-1} sign(eigen_i)*<eigen_i|P>*|eigen_i>  
 * 
 * where  |Q> = Q/Sqrt(Q^2)|S> and thus |Q'>= Q/Sqrt(Q^2)|P> 
 *
 *  
 *  Author: M.Papinutto 
 *  Date: 11.03.2003
 *
 * void sub_lowest_eigenvalues(spinor * const Q, spinor * const P, const int n)
 * 
 * computes: Q=Q-sum_i lambda_i |eigen_i><eigen_i|P>
 *           where eigen_i is the i-th lowest eigenvector and
 *           lambda_i the i-th eigenvalue (of Q^2)
 * Input:
 *   P
 *   n : number of eigenvectors to be subtracted
 * Inout:
 *   Q
 *
 * void assign_add_invert_subtracted_part(spinor * const Q, spinor * const P, const int n)
 *
 * computes: Q = Q + sum_i 1/lambda_i |eigen_i><eigen_i|P>
 *           conventions as obove
 *
 * Input:
 *   P
 *   n : number of eigenvectors to be subtracted
 * Inout:
 *   Q
 *
 * For the last two routines a previous call of
 * eigenvalues or eigenvalues_for_cg must
 * be done
 *
 * Autor: Carsten Urbach <urbach@ifh.de>
 *
 ********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "linalg_eo.h"
#include "eigenvalues.h"
#include "sub_low_ev.h"


/* Q=Q-sum_i lambda_i |eigen_i><eigen_i|P> */
void sub_lowest_eigenvalues(spinor * const Q, spinor * const P, const int n, const int N) {
  int i;
  _Complex double c;
  
  for(i = 0; i < n; i++){
    c = scalar_prod(&(eigenvectors[i*evlength]), P, N, 1);
    c *= -eigenvls[i];
    assign_add_mul(Q, &eigenvectors[i*evlength], c, N);
  }
}

/* Q=P-sum_i |eigen_i><eigen_i|P> */
void assign_sub_lowest_eigenvalues(spinor * const Q, spinor * const P, const int n, const int N) {
  int i;
  _Complex double c;

  assign(Q, P, N);
  
  for(i = 0; i < n; i++){
    c = scalar_prod(&(eigenvectors[i*evlength]), P, N, 1);
    c = -c;
    assign_add_mul(Q, &eigenvectors[i*evlength], c, N);

  }
}

/* Q = Q + sum_i 1/lambda_i |eigen_i><eigen_i|P> */
void assign_add_invert_subtracted_part(spinor * const Q, spinor * const P, const int n, const int N) {
  int i=0;
  _Complex double c;
  double rev=0;

  for(i = 0; i < n; i++){
    c = scalar_prod(&eigenvectors[i*evlength], P, N, 1);
    rev = 1./eigenvls[i];
    c *= rev;
    assign_add_mul(Q, &eigenvectors[i*evlength], c, N);
  }
}

void invert_eigenvalue_part(spinor * const Q, spinor * const P, const int n, const int N) {
  _Complex double c;
  double rev=0;

  assign(Q, P, N);
  for(int i = 0; i < n; ++i)
  {
    c = scalar_prod(&eigenvectors[i*evlength], P, N, 1);
    c *= -inv_eigenvls[i];
    assign_add_mul(Q, &eigenvectors[i*evlength], c, N);
  }
}

