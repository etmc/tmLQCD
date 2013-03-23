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
/*******************************************************************************
 *
 * File gamma.c
 *
 *   void gammaXY ( const Q,  const P)
 *     Makes (*Q) = gammaXY*(*P)   there are 4 gamma_mu, gamma_5 and 4 gamma_5*gamma_mu 
 *
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3spinor.h"
#include "gamma.h"
#ifdef OMP
#include <omp.h>
#endif

/* (*Q) = gammaXY*(*P) */

void gamma0( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < V; ix++){
    _gamma0(g_spinor_field[Q][ix], g_spinor_field[P][ix]);
  }
}
void gamma1( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma1(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma2( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma2(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma3( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma3(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma5(spinor * const l, spinor * const k, const int V){
#ifdef OMP
#pragma omp parallel
  {
#endif
  int ix;
  spinor *r,*s;
#ifdef OMP
#pragma omp for
#endif
  for (ix = 0; ix < V; ix++){
    r=l+ix;
    s=k+ix;
    _vector_assign((*r).s0,(*s).s0);
    _vector_assign((*r).s1,(*s).s1);
    _vector_minus_assign((*r).s2,(*s).s2);
    _vector_minus_assign((*r).s3,(*s).s3);
  }
#ifdef OMP
  } /*OpenMP closing brace */
#endif
}
void gamma5new(spinor * const Q, spinor * const P, const int V){ 
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){ 
    _gamma5(Q[ix], P[ix]); 
  } 
}
void gamma50( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
    for (int ix=0;ix<V;ix++){
    _gamma50(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma51( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma51(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma52( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma52(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}
void gamma53( const int Q,  const int P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix=0;ix<V;ix++){
    _gamma53(g_spinor_field[Q][ix],g_spinor_field[P][ix]);
  }
}

void P_plus(spinor * const Q, spinor * const P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < V; ix++){
    _P_plus(*(Q + ix),*(P + ix));
  }
}

void P_minus(spinor * const Q, spinor * const P, const int V){
#ifdef OMP
#pragma omp parallel for
#endif
  for (int ix = 0; ix < V; ix++){
    _P_minus(*(Q + ix),*(P + ix));
  }
}

void Proj(spinor * const Q, spinor * const P, const int V, const int flag){
  if(flag == 0){ 
#ifdef OMP
#pragma omp parallel for
#endif
    for (int ix = 0; ix < V; ix++){
      _P_plus(*(Q + ix),*(P + ix));
    }
  }
  else if(flag == 1){
#ifdef OMP
#pragma omp parallel for
#endif
    for (int ix = 0; ix < V; ix++){
      _P_minus(*(Q + ix),*(P + ix));
    }
  }
  else{
    printf("wrong flag to Proj! (0/1) \n") ;
  }
}
