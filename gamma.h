/* $Id$  */

#ifndef _GAMMA_H
#define _GAMMA_H

#include "su3.h"

/* Makes (*Q) = gammaXY*(*P)   there are 4 gamma_mu, gamma_5 and 4 gamma_5*gamma_mu  */

void gamma0(const int Q,  const int P, const int V);
void gamma1( const int Q,  const int P, const int V);
void gamma2( const int Q,  const int P, const int V);
void gamma3( const int Q,  const int P, const int V);

/* void gamma5( const int Q,  const int P, const int V); */

void gamma50( const int Q,  const int P, const int V);
void gamma51( const int Q,  const int P, const int V);
void gamma52( const int Q,  const int P, const int V);
void gamma53( const int Q,  const int P, const int V);

#endif
