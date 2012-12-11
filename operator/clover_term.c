/***********************************************************************
 *
 * Copyright (C) 1995 Ulli Wolff, Stefan Sint
 *               2001,2005 Martin Hasenbusch
 *               2011,2012 Carsten Urbach
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
#ifdef SSE
# undef SSE
#endif
#ifdef SSE2
# undef SSE2
#endif
#ifdef SSE3
# undef SSE3
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "operator/clovertm_operators.h"
#include "operator/clover_leaf.h"

// the clover term is written as
//
//   1 + T_{xa\alpha,yb\beta} 
// = 1 + i csw kappa/2 sigma_munu^alphabeta F_munu^ab(x)delta_xy
//
// see hep-lat/9603008 for all glory details
//
// per site we have to store two six-by-six complex matrices.
// As the off-diagonal 3x3 matrices are just inverse to
// each other, we get away with two times three 3x3 complex matrices
//
// these are stored in the array sw[VOLUME][3][2] of type su3
// where x is the space time index
// a runs from 0 to 2
// b runs from 0 to 1
// sw[x][0][0] is the upper diagonal 3x3 matrix 
// sw[x][1][0] the upper off-diagnoal 3x3 matrix
// sw[x][2][0] the lower diagonal 3x3 matrix
// the lower off-diagonal 3x3 matrix would be the inverser of sw[x][1][0]
// 
// identical convention for the second six-by-six matrix
// just with second index set to 1
//
// so the application of the clover term 
// plus twisted mass term to a spinor would just be
// 
// r_0 = sw[0][0] s_0 + sw[1][0] s_1 + i mu s_0
// r_1 = sw[1][0]^-1 s_0 + sw[2][0] s_1 + i mu s_1
// r_2 = sw[0][1] s_2 + sw[1][1] s_3 - i mu s_2
// r_3 = sw[1][1]^-1 s_2 + sw[2][1] s_3 - i mu s_3
//
// suppressing space-time indices

void sw_term(const su3 ** const gf, const double kappa, const double c_sw) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int k,l;
  int x,xpk,xpl,xmk,xml,xpkml,xplmk,xmkml;
  const su3 *w1,*w2,*w3,*w4;
  double ka_csw_8 = kappa*c_sw/8.;
  su3 ALIGN v1,v2,plaq;
  su3 ALIGN fkl[4][4];
  su3 ALIGN magnetic[4],electric[4];
  su3 ALIGN aux;
  

  /*  compute the clover-leave */
  /*  l  __   __
        |  | |  |
        |__| |__|
         __   __
        |  | |  |
        |__| |__| k  */
  
#ifdef OMP
#pragma omp for
#endif
  for(x = 0; x < VOLUME; x++) {
    for(k = 0; k < 4; k++) {
      for(l = k+1; l < 4; l++) {
	xpk=g_iup[x][k];
	xpl=g_iup[x][l];
	xmk=g_idn[x][k];
	xml=g_idn[x][l];
	xpkml=g_idn[xpk][l];
	xplmk=g_idn[xpl][k];
	xmkml=g_idn[xml][k];
	w1=&gf[x][k];
	w2=&gf[xpk][l];
	w3=&gf[xpl][k];
	w4=&gf[x][l];
	_su3_times_su3(v1,*w1,*w2);
	_su3_times_su3(v2,*w4,*w3);
	_su3_times_su3d(plaq,v1,v2);
	w1=&gf[x][l];
	w2=&gf[xplmk][k];
	w3=&gf[xmk][l];
	w4=&gf[xmk][k];
	_su3_times_su3d(v1,*w1,*w2);
	_su3d_times_su3(v2,*w3,*w4);
	_su3_times_su3_acc(plaq,v1,v2);
	w1=&gf[xmk][k];
	w2=&gf[xmkml][l];
	w3=&gf[xmkml][k];
	w4=&gf[xml][l];
	_su3_times_su3(v1,*w2,*w1);
	_su3_times_su3(v2,*w3,*w4);
	_su3d_times_su3_acc(plaq,v1,v2);
	w1=&gf[xml][l];
	w2=&gf[xml][k];
	w3=&gf[xpkml][l];
	w4=&gf[x][k];
	_su3d_times_su3(v1,*w1,*w2);
	_su3_times_su3d(v2,*w3,*w4);
	_su3_times_su3_acc(plaq,v1,v2);
	_su3_dagger(v2,plaq); 
	_su3_minus_su3(fkl[k][l],plaq,v2);
      }
    }

    // this is the one in flavour and colour space
    // twisted mass term is treated in clover, sw_inv and
    // clover_gamma5 and the corresponding nd versions
    _su3_one(sw[x][0][0]);
    _su3_one(sw[x][2][0]);
    _su3_one(sw[x][0][1]);
    _su3_one(sw[x][2][1]);
    
    for(k = 1; k < 4; k++)
    {
      _su3_assign(electric[k], fkl[0][k]);
    }
    _su3_assign(magnetic[1], fkl[2][3]);
    _su3_minus_assign(magnetic[2], fkl[1][3]);
    _su3_assign(magnetic[3], fkl[1][2]);
    
    /*  upper left block 6x6 matrix  */
    
    _itimes_su3_minus_su3(aux,electric[3],magnetic[3]);
    _su3_refac_acc(sw[x][0][0],ka_csw_8,aux);
    
    _itimes_su3_minus_su3(aux,electric[1],magnetic[1]);
    _su3_minus_su3(v2,electric[2],magnetic[2]); 
    _su3_acc(aux,v2);
    _real_times_su3(sw[x][1][0],ka_csw_8,aux);
    
    _itimes_su3_minus_su3(aux,magnetic[3],electric[3]);
    _su3_refac_acc(sw[x][2][0],ka_csw_8,aux);

    /*  lower right block 6x6 matrix */
    
    _itimes_su3_plus_su3(aux,electric[3],magnetic[3]);
    _su3_refac_acc(sw[x][0][1],(-ka_csw_8),aux);

    _itimes_su3_plus_su3(aux,electric[1],magnetic[1]);
    _su3_plus_su3(v2,electric[2],magnetic[2]); 
    _su3_acc(aux,v2);
    _real_times_su3(sw[x][1][1],(-ka_csw_8),aux);

    _itimes_su3_plus_su3(aux,magnetic[3],electric[3]);
    _su3_refac_acc(sw[x][2][1],ka_csw_8,aux);
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}
