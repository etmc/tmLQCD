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

#ifndef _GAMMA_H
#define _GAMMA_H

#include "su3.h"

/* Makes (*Q) = gammaXY*(*P)   there are 4 gamma_mu, gamma_5 and 4 gamma_5*gamma_mu  */

void gamma0(const int Q,  const int P, const int V);
void gamma1( const int Q,  const int P, const int V);
void gamma2( const int Q,  const int P, const int V);
void gamma3( const int Q,  const int P, const int V);

void gamma5(spinor * const Q, spinor * const P, const int V); 

void gamma50( const int Q,  const int P, const int V);
void gamma51( const int Q,  const int P, const int V);
void gamma52( const int Q,  const int P, const int V);
void gamma53( const int Q,  const int P, const int V);

void P_plus(spinor * const Q, spinor * const P, const int V); 
void P_minus(spinor * const Q, spinor * const P, const int V); 
void Proj(spinor * const Q, spinor * const P, const int V, const int flag); 

#endif



