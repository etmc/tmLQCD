/*****************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Deflating CG using eigenvectors computed using ARPACK
 *
 * Author: A.M. Abdel-Rehim (amabdelrehim@gmail.com)
 *         November 2014
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
 *
 ****************************************************************************/
/* A sample input is given in the sample-input folder */

#ifndef _DEFLATOR_H
#define _DEFLATOR_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/eigenvalues_arpack.h"


typedef struct deflator_params_t {

  char type_name[100];
  int type;

  int eoprec;

  matrix_mult f;    /* f(s,r)   computes s=A*r, i.e. matrix-vector multiply in double precision */
  matrix_mult32 f32;  /* f32(s,r) computes s=A*r, i.e. matrix-vector multiply in single precision */

  matrix_mult f_final;    /* (IN) final operator application during projection of type 1 */
  matrix_mult f_initial;  /* (IN) initial operator application during projection of type 1 */

  int projection_type;

  void *evecs;
  void *evals;
  int prec;

  /**********************************
   * arpack parameters
   **********************************/

   int    nconv;         /* number of converged eigenvectors */
   int    nev;           /* number of eigenvectors to be computed by arpack*/
   int    ncv;           /* Size of the subspace used by arpack with the condition (nev+1) =< ncv */
   int    evals_kind;    /* type of eigenvalues to be computed 
                                     0 eigenvalues of smallest real part
                                     1 eigenvalues of largest real part 
                                     2 eigenvalues of smallest absolute value
                                     3 eigenvalues of largest absolute value*/
   int    comp_evecs;    /* 0 don't compute the resiudals of the eigenvalues
                                     1 compute the residulas of the eigenvalues*/
   double eig_tol;       /* tolerance for computing eigenvalues with arpack */
   int    eig_maxiter;   /* maximum number of iterations to be used by arpack*/
   char   logfile[500];  /* file name for the logfile used by arpack*/

  /*************************************
   * Chebychev polynomial acceleraton
   * paprameters
   *************************************/
   int use_acc;      /* type of acceleration to be used:  */
                     /*   0 no acceleration, 1 compute eiegnvectors of the acceleration polynomial,   */
   int cheb_k;       /* order of the polynomial used is k+1 and the lowest value is k=-1 which correspond to T_0 */
   double op_evmin;  /* lowest boundary of the interval for the polynomial acceleration */
   double op_evmax;  /* highest boundary for the interval for the polynomial acceleration */
   
  /*************************************
   * parameters to control of I/O for
   * eigenvectors
   *************************************/
   int write_ev;
   int read_ev;
   char evecs_filename[500];    /* file name for the evecs used by arpack*/
   char evecs_fileformat[200];  /* file format for evecs used by arpack */

   int evecs_writeprec;

  /*************************************
   * init function
   *************************************/
  int (*init)(struct deflator_params_t*);
  int (*fini)(struct deflator_params_t*);

} deflator_params_t;

int make_exactdeflator( deflator_params_t *deflator_params);

#endif
