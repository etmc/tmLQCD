/***********************************************************************
 * $Id$ 
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "global.h"
#include "default_input_values.h"
#include "su3.h"
#include "operator.h"

void dummy_D(spinor * const, spinor * const);
int dummy_inverter(spinor * const, spinor * const, spinor * const, spinor * const,
		   const int op_id);

operator operator_list[max_no_operators];

int no_operators = 0;

int add_operator(const int type) {

  if(no_operators == max_no_operators) {
    fprintf(stderr, "maximal number of operators %d exceeded!\n", max_no_operators);
    exit(-1);
  }
  operator_list[no_operators].type = type;
  operator_list[no_operators].kappa = _default_g_kappa;
  operator_list[no_operators].mu = _default_g_mu;
  operator_list[no_operators].s = 0.;
  operator_list[no_operators].m = 0.;
  operator_list[no_operators].sloppy_precision = _default_g_sloppy_precision_flag;
  operator_list[no_operators].coefs = NULL;
  if(operator_list[no_operators].type == OVERLAP) {
    operator_list[no_operators].even_odd = 0;
    operator_list[no_operators].solver = 13;
  }
  else {
    operator_list[no_operators].even_odd = _default_even_odd_flag;
    operator_list[no_operators].solver = _default_solver_flag;
  }
  operator_list[no_operators].applyD = &dummy_D;
  operator_list[no_operators].applyQ = &dummy_D;
  operator_list[no_operators].applyQp = &dummy_D;
  operator_list[no_operators].applyQm = &dummy_D;
  operator_list[no_operators].applyQsq = &dummy_D;
  operator_list[no_operators].inverter = &dummy_inverter;

  operator_list[no_operators].initialised = 1;

  no_operators++;
  return(no_operators);
}

int init_operators() {
  int i;
  for(i = 0; i < no_operators; i++) {
    /* do something here for all operators */
  }  
  return(0);
}

void dummy_D(spinor * const s, spinor * const r) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_D was called. Was that really intended?\n");
  } 
  return;
}

int dummy_inverter(spinor * const p, spinor * const q, spinor * const r, spinor * const s,
		   const int op_id) {
  if(g_proc_id == 0) {
    fprintf(stderr, "dummy_inverter was called. Was that really intended?\n");
  } 
  return(0);
}
