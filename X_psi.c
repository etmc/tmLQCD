/***********************************************************************
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
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "operator/D_psi.h"
#include "gamma.h"
#include "X_psi.h"
#include "operator/tm_operators.h"
#include "solver/solver.h"
#include "read_input.h"

void DdaggerD_plus_M(spinor * const R, spinor * const S)
{
  double g_muWithoutMStarSquare=g_mu;
  g_mu=sqrt(g_mu*g_mu+mstarsq);
  Q_pm_psi(R, S);
  g_mu=g_muWithoutMStarSquare;
  
/*  spinor *aux_ = NULL, *aux;
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

  // we have to apply DdagerD and M*^2 to the same field S
  twokmu=g_mu;
  g_mu=0.;
  //org: 
  //D_psi(R, S);
  //gamma5(aux, R, VOLUME);
  //D_psi(R, aux);
  //gamma5(R, R, VOLUME);
  D_psi(aux, S);
  gamma5(R, aux, VOLUME);
  D_psi(aux,R);
  gamma5(R, aux, VOLUME);
  
  g_mu=twokmu;
  g_musq=g_mu*g_mu;
  assign_add_mul_r(R, aux2, mstarsq, VOLUME);
  if(g_musq!=0) assign_add_mul_r(R, aux2, g_musq, VOLUME);
  
  free(aux_);
  free(aux2_);*/
}

#define X_psiSIterations 5000
#define X_psiSPrecision 1.e-6


void X_psi(spinor * const R, spinor * const S, double const mstarsq){

  //  double a = -2*mstar*mstar;
  double a = -2*mstarsq;
  double b = 1.;
  double g_muWithoutMStarSquare=g_mu;

  /*cg_her(out spinor, in spinor, max iter, solver precision, flag relative precision default 0, volume, operator to invert)*/
  if(g_proc_id == 0) printf("Using CPU for inversion\n");
  cg_her( R, S, X_psiSIterations, X_psiSPrecision, 0, VOLUME, &DdaggerD_plus_M);

/*//// Test
  spinor *aux_ = NULL, *aux;
  int N = VOLUMEPLUSRAND;

#if ( defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif

  //Q_pm_psi_gpu(aux,R);
  DdaggerD_plus_M(aux,R);
  diff(aux,S,aux,N);
  double t=square_norm(aux,N,1);
  printf("TestMStar %lf\n",t);
  exit(1);
*//// Test
  assign_mul_add_mul_r( R, S, a, b, VOLUME);
}


void X_psiSquare(spinor * const R, spinor * const S, double const mstarsq)
{//inverts DD^+DD^+ instead of DD^+ but performs poorly
  spinor *aux_,*aux;
  {
    #if ( defined SSE || defined SSE2 || defined SSE3 )
      aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
      aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
    #else
      aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
      aux = aux_;
    #endif
  }

  X_psi(aux, S, mstarsq);
  X_psi(R, aux, mstarsq);

  free(aux_);
}
 
