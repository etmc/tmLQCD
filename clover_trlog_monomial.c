/***********************************************************************
 *
 * Copyright (C) 1008 Carsten Urbach
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
#include <time.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "su3spinor.h"
#include "ranlxd.h"
#include "start.h"
#include "read_input.h"
#include "clover_leaf.h"

#include "monomial.h"
#include "clover_trlog_monomial.h"

void clover_trlog_derivative(const int no, hamiltonian_field_t * const hf) {

  return;
}


void clover_trlog_heatbath(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  /*compute the contribution from the clover-term*/
  mnl->energy0 = 2.*sw_trace(1);
  return;
}

double clover_trlog_acc(const int id, hamiltonian_field_t * const hf) {
  monomial * mnl = &monomial_list[id];
  sw_term(hf->gaugefield, mnl->kappa, mnl->c_sw); 
  /*compute the contribution from the clover-term*/
  mnl->energy1 = 2.*sw_trace(1);  
  return(mnl->energy1 - mnl->energy0);
}
