/***********************************************************************
 * $Id: det_monomial.c 1764 2011-04-21 14:16:16Z deuzeman $
 *
 * Copyright (C) 2011 Elena Garcia-Ramos
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
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "D_psi.h"
#include "gamma.h"
#include "X_psi.h"
#include "tm_operators.h"
#include "read_input.h"

void DdaggerD_plus_M(spinor * const R, spinor * const S)
{
  spinor *aux_ = NULL, *aux;
  spinor *aux2_ = NULL, *aux2;
  int N = VOLUMEPLUSRAND;
  double twokmu, g_musq;

#if ( defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
  aux2_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux2 = (spinor *)(((unsigned long int)(aux2_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
  aux2_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux2 = aux2_;
#endif

  assign(aux2, S, VOLUME);

  /* we have to apply DdagerD and M*^2 to the same field S*/
  twokmu=g_mu;
  g_mu=0.;
  D_psi(R, S);
  gamma5(aux, R, VOLUME);
  D_psi(R, aux);
  gamma5(R, R, VOLUME);
  
  g_mu=twokmu;
  g_musq=g_mu*g_mu;
  assign_add_mul_r(R, aux2, mstarsq, VOLUME);
  assign_add_mul_r(R, aux2, g_musq, VOLUME);
  
  free(aux_);
  free(aux2_);
}




void X_psi(spinor * const R, spinor * const S, double const mstarsq){

  //  double a = -2*mstar*mstar;
  double a = -2*mstarsq;
  double b = 1.;

  /*cg_her(out spinor, in spinor, max iter, solver precision, flag relative precision default 0, volume, operator to invert)*/
  cg_her( R, S, 5000, 1.e-12, 0, VOLUME, &DdaggerD_plus_M);
  
  assign_mul_add_mul_r( R, S, a, b, VOLUME);
  
}
 
