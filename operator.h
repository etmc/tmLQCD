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
#include "operator_types.h"
#include "misc_types.h"

#define max_no_operators 10

typedef struct {
  /* ID of the operator */
  int type;
  int id;
  /* for overlap */
  int n_cheby;
  int deg_poly;
  int no_ev;
  
  SloppyPrecision sloppy_precision;
  int even_odd_flag;
  int solver;
  int N_s;
  int initialised;
  int rel_prec;
  int maxiter;
  int iterations;
  int prop_precision;
  int write_prop_flag;
  int no_flavours;
  int DownProp;
  int no_ev_index;
  ExternalInverter external_inverter;
  CompressionType compression_type;

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
  /* ---------------------------------------*/
  //even-even part of the even-odd operator 
  void (*applyMee) (spinor * const, spinor * const, double const);
  //inverse of the even-even part of the even-odd operator 
  void (*applyMeeInv) (spinor * const, spinor * const, double const);
  //full operator M acting on a spinor given as even and odd parts separately
  void (*applyM) (spinor * const, spinor * const, spinor * const, spinor * const); 
  //full operator Q=gamma5*M on a spinor given as even and odd parts separately
  void (*applyQ) (spinor * const, spinor * const, spinor * const, spinor * const); 
  //either: the full operator Q^+ on lexiographic spinor 
  //or    : eo-preconditioned Q^+ on odd part of an eo ordered spinor
  void (*applyQp) (spinor * const, spinor * const); 
  //either : the full operator Q^- on lexiographic spinor 
  //or     : eo-preconditioned Q^- on odd part of an eo ordered spinor
  void (*applyQm) (spinor * const, spinor * const); 
  //either: the full operator Q^+*Q^- on lexiographic spinor 
  //or    : eo-preconditioned Q^+*Q^- on odd part of an eo ordered spinor
  void (*applyQsq) (spinor * const, spinor * const);
  //either: the full operator M^+ on lexiographic spinor 
  //or    : eo-preconditioned M^+ on odd part of an eo ordered spinor
  void (*applyMp) (spinor * const, spinor * const); 
  //either: the full operator M^- on lexiographic spinor 
  //or    : eo-preconditioned M^- on odd part of an eo ordered spinor
  void (*applyMm) (spinor * const, spinor * const); 
  //EO preconditoned Hermitian operator for the non-degenerate doublet (more explanantion needed here).
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

void op_set_globals(const int op_id);
void op_backup_restore_globals(const backup_restore_t mode);

#endif
