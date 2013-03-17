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

// now we sum up all term from the clover term
// after sw_spinor and sw_deriv have been called

void sw_all(hamiltonian_field_t * const hf, const double kappa, 
	    const double c_sw) {
#ifdef OMP
#pragma omp parallel
  {
#endif

  int k,l;
  int x,xpk,xpl,xmk,xml,xpkml,xplmk,xmkml;
  const su3 *w1,*w2,*w3,*w4;
  double ka_csw_8 = kappa*c_sw/8.;
  su3 ALIGN v1,v2,vv1,vv2,plaq;
  su3 ALIGN vis[4][4];

#ifdef OMP
#pragma omp for
#endif
  for(x = 0; x < VOLUME; x++) {
    _minus_itimes_su3_plus_su3(vis[0][1],swm[x][1],swm[x][3]);
    _su3_minus_su3(vis[0][2],swm[x][1],swm[x][3]);
    _itimes_su3_minus_su3(vis[0][3],swm[x][2],swm[x][0]);
    
    _minus_itimes_su3_plus_su3(vis[2][3],swp[x][1],swp[x][3]);
    _su3_minus_su3(vis[1][3],swp[x][3],swp[x][1]);
    _itimes_su3_minus_su3(vis[1][2],swp[x][2],swp[x][0]);

    // project to the traceless anti-hermitian part
    _su3_dagger(v1,vis[0][1]); 
    _su3_minus_su3(vis[0][1],vis[0][1],v1);
    _su3_dagger(v1,vis[0][2]); 
    _su3_minus_su3(vis[0][2],vis[0][2],v1);
    _su3_dagger(v1,vis[0][3]); 
    _su3_minus_su3(vis[0][3],vis[0][3],v1);
    _su3_dagger(v1,vis[2][3]); 
    _su3_minus_su3(vis[2][3],vis[2][3],v1);
    _su3_dagger(v1,vis[1][3]); 
    _su3_minus_su3(vis[1][3],vis[1][3],v1);
    _su3_dagger(v1,vis[1][2]); 
    _su3_minus_su3(vis[1][2],vis[1][2],v1);
    
    for(k = 0; k < 4; k++) {
      for(l = k+1; l < 4; l++) {
	xpk=g_iup[x][k];
	xpl=g_iup[x][l];
	xmk=g_idn[x][k];
	xml=g_idn[x][l];
	xpkml=g_idn[xpk][l];
	xplmk=g_idn[xpl][k];
	xmkml=g_idn[xml][k];
	w1=&hf->gaugefield[x][k];
	w2=&hf->gaugefield[xpk][l];
	w3=&hf->gaugefield[xpl][k];   /*dag*/
	w4=&hf->gaugefield[x][l];     /*dag*/
	
	_su3_times_su3(v1,*w1,*w2);
	_su3_times_su3(v2,*w4,*w3);
	_su3_times_su3d(plaq,v1,v2);
	
	_su3_times_su3(vv1,plaq,vis[k][l]);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[x][k], -2.*ka_csw_8, vv1);

	_su3d_times_su3(vv2,*w1,vv1); 
	_su3_times_su3(vv1,vv2,*w1);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xpk][l], -2.*ka_csw_8, vv1);
	
	_su3_times_su3(vv2,vis[k][l],plaq); 
	_su3_dagger(vv1,vv2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[x][l], -2.*ka_csw_8, vv1);

	_su3d_times_su3(vv2,*w4,vv1); 
	_su3_times_su3(vv1,vv2,*w4);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xpl][k], -2.*ka_csw_8, vv1);
	
	w1=&hf->gaugefield[x][l];
	w2=&hf->gaugefield[xplmk][k];   /*dag*/
	w3=&hf->gaugefield[xmk][l];     /*dag*/
	w4=&hf->gaugefield[xmk][k];
	_su3_times_su3d(v1,*w1,*w2);
	_su3d_times_su3(v2,*w3,*w4);
	_su3_times_su3(plaq,v1,v2);
	
	_su3_times_su3(vv1,plaq,vis[k][l]);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[x][l], -2.*ka_csw_8, vv1);
	
	_su3_dagger(vv1,v1); 
	_su3_times_su3d(vv2,vv1,vis[k][l]);
	_su3_times_su3d(vv1,vv2,v2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xplmk][k], -2.*ka_csw_8, vv1);

	_su3_times_su3(vv2,*w3,vv1); 
	_su3_times_su3d(vv1,vv2,*w3);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xmk][l], -2.*ka_csw_8, vv1);

	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xmk][k], -2.*ka_csw_8, vv2);
	
	w1=&hf->gaugefield[xmk][k];   /*dag*/
	w2=&hf->gaugefield[xmkml][l]; /*dag*/
	w3=&hf->gaugefield[xmkml][k];
	w4=&hf->gaugefield[xml][l];
	_su3_times_su3(v1,*w2,*w1);
	_su3_times_su3(v2,*w3,*w4);
	
	_su3_times_su3d(vv1,*w1,vis[k][l]);
	_su3_times_su3d(vv2,vv1,v2);
	_su3_times_su3(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xmk][k], -2.*ka_csw_8, vv1);

	_su3_times_su3(vv2,*w2,vv1); 
	_su3_times_su3d(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xmkml][l], -2.*ka_csw_8, vv1);

	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xmkml][k], -2.*ka_csw_8, vv2);

	_su3d_times_su3(vv1,*w3,vv2); 
	_su3_times_su3(vv2,vv1,*w3);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xml][l], -2.*ka_csw_8, vv2);
	
	w1=&hf->gaugefield[xml][l];   /*dag*/
	w2=&hf->gaugefield[xml][k];
	w3=&hf->gaugefield[xpkml][l];
	w4=&hf->gaugefield[x][k];     /*dag*/
	_su3d_times_su3(v1,*w1,*w2);
	_su3_times_su3d(v2,*w3,*w4);
	
	_su3_times_su3d(vv1,*w1,vis[k][l]);
	_su3_times_su3d(vv2,vv1,v2);
	_su3_times_su3d(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xml][l], -2.*ka_csw_8, vv1);
	
	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xml][k], -2.*ka_csw_8, vv2);

	_su3d_times_su3(vv1,*w2,vv2); 
	_su3_times_su3(vv2,vv1,*w2);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[xpkml][l], -2.*ka_csw_8, vv2);

	_su3_dagger(vv2,v2);  
	_su3_times_su3d(vv1,vv2,v1);
	_su3_times_su3d(vv2,vv1,vis[k][l]);
 	_trace_lambda_mul_add_assign_nonlocal(hf->derivative[x][k], -2.*ka_csw_8, vv2);
      }
    }
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}
