/*******************************************************************************
 * $Id$
 *
 * File gamma.c
 *
 *   void gammaXY ( const Q,  const P)
 *     Makes (*Q) = gammaXY*(*P)   there are 4 gamma_mu, gamma_5 and 4 gamma_5*gamma_mu 
 *
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "global.h"
#include "su3spinor.h"
#include "gamma.h"

/* (*Q) = gammaXY*(*P) */

void gamma0( const int Q,  const int P, const int V){
  int ix;
  
  for (ix = 0; ix < V; ix++){
    _gamma0(spinor_field[Q][ix], spinor_field[P][ix]);
  }
}
void gamma1( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma1(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
void gamma2( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma2(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
void gamma3( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma3(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
/* void gamma5( const int Q,  const int P, const int V){ */
/*   int ix; */
  
/*   for (ix=0;ix<V;ix++){ */
/*     _gamma5(spinor_field[Q][ix],spinor_field[P][ix]); */
/*   } */
/* } */
void gamma50( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma50(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
void gamma51( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma51(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
void gamma52( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma52(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
void gamma53( const int Q,  const int P, const int V){
  int ix;
  
  for (ix=0;ix<V;ix++){
    _gamma53(spinor_field[Q][ix],spinor_field[P][ix]);
  }
}
