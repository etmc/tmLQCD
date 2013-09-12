
/***********************************************************************
 *
 * Copyright (C) 2009 Carsten Urbach
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

#ifndef _OPERATOR_H
#define _OPERATOR_H

#include <io/utils.h>
#include "solver/dirac_operator_eigenvectors.h"
#include "su3.h"
#include "solver/solver_params.h"
                                    

#define TMWILSON 0
#define OVERLAP 1
#define WILSON 2
#define DBTMWILSON 3
#define CLOVER 4
#define DBCLOVER 5

#define max_no_operators 10

typedef struct {
  /* ID of the operator */
  int type;
  int id;
  /* for overlap */
  int n_cheby;
  int deg_poly;
  int no_ev;
  
  int sloppy_precision;
  int even_odd_flag;
  int solver;
  int N_s;
  int initialised;
  int rel_prec;
  int maxiter;
  int iterations;
  int prop_precision;
  int no_flavours;
  int DownProp;
  int no_ev_index;

  int error_code;

  double kappa;
  /* for twisted */
  double mu;
  double mubar;
  /* for 2 flavour twisted */
  double epsbar;
  /* solver residue */
  double eps_sq;
  /* clover coefficient */
  double c_sw;
  /* precision reached during inversion */
  double reached_prec;
  /* for the overlap */
  double m;
  double s;
  double ev_qnorm;
  double ev_minev;
  double ev_prec;
  int ev_readwrite;
  /* generic place for sources */
  spinor *sr0, *sr1, *sr2, *sr3;
  /* generic place for propagators */
  spinor *prop0, *prop1, *prop2, *prop3;

  /*solver parameters struct*/
  solver_params_t solver_params;

  /* multiple masses for CGMMS */
  double extra_masses[MAX_EXTRA_MASSES];
  int no_extra_masses;


  /* chebyshef coefficients for the overlap */
  double * coefs;
  /* various versions of the Dirac operator */
  void (*applyM) (spinor * const, spinor * const);
  void (*applyQ) (spinor * const, spinor * const);
  /* with even/odd */
  void (*applyQp) (spinor * const, spinor * const);
  void (*applyQm) (spinor * const, spinor * const);
  void (*applyQsq) (spinor * const, spinor * const);
  void (*applyMp) (spinor * const, spinor * const);
  void (*applyMm) (spinor * const, spinor * const);
  void (*applyDbQsq) (spinor * const, spinor * const, spinor * const, spinor * const);
  /* the generic invert function */
  void (*inverter) (const int op_id, const int index_start, const int write_prop);
  /* write the propagator */
  void (*write_prop) (const int op_id, const int index_start, const int append_);
  char * conf_input;

  spinorPrecWS *precWS;
  
} operator;

/* operator list defined in operator.c */
extern operator operator_list[max_no_operators];
extern int no_operators;

int add_operator(const int type);
int init_operators();

#endif
