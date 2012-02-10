/***********************************************************************
 *
 * Copyright (C) 1995 Ulli Wolff, Stefan Sint
 *               2001,2005 Martin Hasenbusch
 *               2011 Carsten Urbach
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
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "clover.h"
#include "clover_leaf.h"

const double tiny_t = 1.0e-20;

su3 ** swm, ** swp;

void sw_term(su3 ** const gf, const double kappa, const double c_sw) {
  int k,l;
  int x,xpk,xpl,xmk,xml,xpkml,xplmk,xmkml;
  su3 *w1,*w2,*w3,*w4;
  double ka_csw_8 = kappa*c_sw/8.;
  static su3 v1,v2,plaq;
  static su3 fkl[4][4];
  static su3 magnetic[4],electric[4];
  static su3 aux;
  
  /*  compute the clover-leave */
  /*  l  __   __
        |  | |  |
        |__| |__|
         __   __
        |  | |  |
        |__| |__| k  */
  
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
    
    /*  upper left block matrix  */
    
    _itimes_su3_minus_su3(aux,electric[3],magnetic[3]);
    _su3_refac_acc(sw[x][0][0],ka_csw_8,aux);
    
    _itimes_su3_minus_su3(aux,electric[1],magnetic[1]);
    _su3_minus_su3(v2,electric[2],magnetic[2]); 
    _su3_acc(aux,v2);
    _real_times_su3(sw[x][1][0],ka_csw_8,aux);
    
    _itimes_su3_minus_su3(aux,magnetic[3],electric[3]);
    _su3_refac_acc(sw[x][2][0],ka_csw_8,aux);

    /*  lower right block matrix */
    
    _itimes_su3_plus_su3(aux,electric[3],magnetic[3]);
    _su3_refac_acc(sw[x][0][1],(-ka_csw_8),aux);

    _itimes_su3_plus_su3(aux,electric[1],magnetic[1]);
    _su3_plus_su3(v2,electric[2],magnetic[2]); 
    _su3_acc(aux,v2);
    _real_times_su3(sw[x][1][1],(-ka_csw_8),aux);

    _itimes_su3_plus_su3(aux,magnetic[3],electric[3]);
    _su3_refac_acc(sw[x][2][1],ka_csw_8,aux);
  }
  return;
}

/*
  !--------------------------------------------------------------!
  !  The subroutine sw_invert is needed for the                  !
  !  even_odd preconditioned Dirac operator with SW improvement. !
  !  Details can be found in  the notes sw.ps on tsun.desy.de    !
  !  by P. Weisz and U. Wolff.                                   !
  !--------------------------------------------------------------!
  !  inversion in place of _Complex double matrix a without pivoting     !
  !  triangularization by householder reflexions                 !
  !  inversion of triangular matrix                              !
  !  inverse reflexions                                          !
  !--------------------------------------------------------------!
  !  a square matrix, dimensioned 0:n-1                          !
  !  itrouble is counted up, when a dangerously small diagonal   !
  !  element is encountered in the tringular matrix              !
  !  has to be initialized outside                               !
  !                                                              !
  !  Author: U. Wolff, adapted to fortran90 by S. Sint, 29/10/95 !
  !--------------------------------------------------------------!
  !  ported to C by M.Hasenbusch Wed Oct 24 15:46:46 MEST 2001   !
  !______________________________________________________________!
*/
#define nm1 5
int six_invert(_Complex double a[6][6])
{
  static _Complex double d[nm1+1],u[nm1+1];
  static _Complex double sigma,z;
  static double p[nm1+1];
  static double s,q;
  int i,j,k;
  int ifail;
  ifail=0;
  for(k = 0; k < nm1; ++k)
  {
    s=0.0;
    for(j = k+1; j <= nm1; ++j)
      s += conj(a[j][k]) * a[j][k];
    s = sqrt(1. + s / (conj(a[k][k]) * a[k][k]));
    sigma = s * a[k][k];

    a[k][k] += sigma;
    p[k] = creal(conj(sigma) * a[k][k]);
    q = conj(sigma) * sigma;
    if (q < tiny_t)
      ifail++;
    d[k] = -conj(sigma) / q;

    /* reflect all columns to the right */
    for(j = k+1; j <= nm1; j++)
    {
      z = 0.0;
      for(i = k; i <= nm1; i++)
	z += conj(a[i][k]) * a[i][j];
      z /= p[k];
      for(i = k; i <= nm1; i++)
	a[i][j] -= z * a[i][k];
    }
  }
  sigma = a[nm1][nm1];
  q = conj(sigma) * sigma;
  if (q < tiny_t)
    ifail++;
  d[nm1] = -conj(sigma) / q;

  /*  inversion of upper triangular matrix in place
      (diagonal elements done already): */

  for(k = nm1; k >= 0; k--) {
    for(i = k-1; i >= 0;i--) {
      z = 0.0;
      for(j = i+1; j < k; j++)
	z += a[i][j] * a[j][k];
      z += a[i][k] * d[k];
      a[i][k] = -z * d[i];
    }
  }     
  /* execute reflexions in reverse order from the right: */
  
  a[nm1][nm1] = d[nm1];
  for(k = nm1-1; k >= 0; k--)
  {
    for(j=k;j<=nm1;j++)
      u[j] = a[j][k];
    a[k][k] = d[k];
    for(j = k+1; j <= nm1; j++)
      a[j][k] = 0.0;
    for(i = 0; i <= nm1; i++)
    {
      z = 0.0;
      for(j = k; j <= nm1; j++)
        z += a[i][j] * u[j];
      z /= p[k];         /* normalization */
      
      for(j = k; j <= nm1; j++)
        a[i][j] -= conj(u[j]) * z; /* reflexion */
    }
  }
  return ifail;
}
    
double six_det(_Complex double a[6][6])
{
  static _Complex double sigma,z;
  static _Complex double det;
  static double p[nm1+1];
  static double s,q;
  int i,j,k;
  int ifail;
  ifail=0;
  /* compute the determinant:*/
  det = 1.0;

  for(k = 0; k < nm1; k++)
  {
    s=0.0;
    for(j = k+1; j <= nm1; ++j)
      s += conj(a[j][k]) * a[j][k];
    s = sqrt(1. + s / (conj(a[k][k]) * a[k][k]));
    sigma = s * a[k][k];

    /* determinant */
    det = sigma * det;
    q   = sigma * conj(sigma);
    if (q < tiny_t)
      ifail++;

    a[k][k] += sigma;
    p[k]     = sigma * conj(a[k][k]);

    /* reflect all columns to the right */
    for(j = k+1; j <= nm1; j++)
    {
      z = 0.;
      for(i = k; i <= nm1; i++)
	z += conj(a[i][k]) * a[i][j];
      z /= p[k];
      for(i = k; i <= nm1; i++)
	a[i][j] -= z * a[i][k];
    }
  }
  sigma = a[nm1][nm1];

  /* determinant */
  det = sigma * det;
  q = conj(sigma) * sigma;

  if(q < tiny_t)
    ifail++;

  return det;
}

/*definitions needed for the functions sw_trace(int ieo) and sw_trace(int ieo)*/
inline void populate_6x6_matrix(_Complex double a[6][6], su3 * C, const int row, const int col) {
  a[0+row][0+col] = C->c00;
  a[0+row][1+col] = C->c01;
  a[0+row][2+col] = C->c02;
  a[1+row][0+col] = C->c10;
  a[1+row][1+col] = C->c11;
  a[1+row][2+col] = C->c12;
  a[2+row][0+col] = C->c20;
  a[2+row][1+col] = C->c21;
  a[2+row][2+col] = C->c22;
  return;
}

inline void get_3x3_block_matrix(su3 * C, _Complex double a[6][6], const int row, const int col) {
  C->c00 = a[0+row][0+col];
  C->c01 = a[0+row][1+col];
  C->c02 = a[0+row][2+col];
  C->c10 = a[1+row][0+col];
  C->c11 = a[1+row][1+col];
  C->c12 = a[1+row][2+col];
  C->c20 = a[2+row][0+col];
  C->c21 = a[2+row][1+col];
  C->c22 = a[2+row][2+col];
  return;
}

/*definitions needed for the functions sw_trace(int ieo) and sw_trace(int ieo)*/
#define _a_C(A, B, C)\
  a[0+(A)][0+(B)]=(C).c00;			\
  a[0+(A)][1+(B)]=(C).c01;			\
  a[0+(A)][2+(B)]=(C).c02;			\
  a[1+(A)][0+(B)]=(C).c10;			\
  a[1+(A)][1+(B)]=(C).c11;			\
  a[1+(A)][2+(B)]=(C).c12;			\
  a[2+(A)][0+(B)]=(C).c20;			\
  a[2+(A)][1+(B)]=(C).c21;			\
  a[2+(A)][2+(B)]=(C).c22;

#define _C_a(A, B, C)\
  (C).c00=a[0+(A)][0+(B)];			\
  (C).c01=a[0+(A)][1+(B)];			\
  (C).c02=a[0+(A)][2+(B)];			\
  (C).c10=a[1+(A)][0+(B)];			\
  (C).c11=a[1+(A)][1+(B)];			\
  (C).c12=a[1+(A)][2+(B)];			\
  (C).c20=a[2+(A)][0+(B)];			\
  (C).c21=a[2+(A)][1+(B)];			\
  (C).c22=a[2+(A)][2+(B)];

double sw_trace(const int ieo) {
  int i,x,icx,ioff;
  static su3 v;
  static _Complex double a[6][6];
  static double tra;
  static double ks,kc,tr,ts,tt;
  
  ks=0.0;
  kc=0.0;
  
  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    for(i=0;i<2;i++) {
      populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
      populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
      _su3_dagger(v, sw[x][1][i]); 
      populate_6x6_matrix(a, &v, 3, 0);
      populate_6x6_matrix(a, &sw[x][2][i], 3, 3);
      tra = log(six_det(a));
      
      tr=tra+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  }
  kc=ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return(ks);
#else
  return(kc);
#endif

}



void sw_invert(const int ieo) {
  int icx,ioff, err=0;
  int i,x;
  static su3 v;
  static _Complex double a[6][6];
  if(ieo==0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }

  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];

    for(i = 0; i < 2; i++) {
      populate_6x6_matrix(a, &sw[x][0][i], 0, 0);
      populate_6x6_matrix(a, &sw[x][1][i], 0, 3);
      _su3_dagger(v, sw[x][1][i]); 
      populate_6x6_matrix(a, &v, 3, 0);
      populate_6x6_matrix(a, &sw[x][2][i], 3, 3);

      err = six_invert(a); 
      /* here we need to catch the error! */
      if(err > 0 && g_proc_id == 0) {
	printf("# %d\n", err);
	err = 0;
      }

      /*  copy "a" back to sw_inv */
      get_3x3_block_matrix(&sw_inv[x][0][i], a, 0, 0);
      get_3x3_block_matrix(&sw_inv[x][1][i], a, 0, 3);
      get_3x3_block_matrix(&sw_inv[x][2][i], a, 3, 3);
    }
  }
  return;
}

void sw_deriv(const int ieo) {
  int ioff,icx;
  int x;
  static su3 lswp[4],lswm[4];
  
  /*  compute the clover-leave */
  /*  l  __   __
        |  | |  |
        |__| |__|
         __   __
        |  | |  |
        |__| |__| k  */
  
  
  /* convention: Tr colver-leaf times insertion */
  /* task : put the matrix in question to the front */
  if(ieo == 0) {
    ioff=0;
  } 
  else {
    ioff = (VOLUME+RAND)/2;
  }
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0],sw_inv[x][0][1],sw_inv[x][0][0]);
    _su3_plus_su3(lswp[1],sw_inv[x][1][1],sw_inv[x][1][0]);
    _su3_plus_su3(lswp[2],sw_inv[x][2][1],sw_inv[x][2][0]);
    _su3_dagger(lswp[3],lswp[1]);
    _su3_assign(lswp[1], lswp[3]);
    _su3_minus_su3(lswm[0],sw_inv[x][0][1],sw_inv[x][0][0]);
    _su3_minus_su3(lswm[1],sw_inv[x][1][1],sw_inv[x][1][0]);
    _su3_minus_su3(lswm[2],sw_inv[x][2][1],sw_inv[x][2][0]);
    _su3_dagger(lswm[3],lswm[1]);
    _su3_assign(lswm[1], lswm[3]);
    
    /* add up the swm[0]  and swp[2] */
    _su3_acc(swm[x][0],lswm[0]);
    _su3_acc(swm[x][1],lswm[1]);
    _su3_acc(swm[x][2],lswm[2]);
    _su3_acc(swm[x][3],lswm[3]);
    _su3_acc(swp[x][0],lswp[0]);
    _su3_acc(swp[x][1],lswp[1]);
    _su3_acc(swp[x][2],lswp[2]);
    _su3_acc(swp[x][3],lswp[3]);
  }
  return;
}


/* ieo : even =0 ; odd =1  k: spinor-field on the right; l: on the left */
void sw_spinor(const int ieo, spinor * const kk, spinor * const ll) {
  int ioff;
  int icx;
  int x;
  spinor *r,*s;
  static su3 v1,v2,v3;
  static su3 u1,u2,u3;
  static su3 lswp[4],lswm[4];
  
  /*  compute the clover-leave */
  /*  l  __   __
        |  | |  |
        |__| |__|
 	 __   __
        |  | |  |
        |__| |__| k  */
  
  
  /* convention: Tr colver-leaf times insertion */
  /* task : put the matrix in question to the front */

  if(ieo == 0) {
    ioff=0;
  } 
  else {
    ioff=(VOLUME+RAND)/2;
  }
  /************************ loop over half of the lattice sites ***********/
  
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++) {
    x = g_eo2lexic[icx];
    r = kk + icx - ioff;
    s = ll + icx - ioff;
    
    /*
      s1= sw(1,1,l,t)     *phi(1,l,t) +sw(2,1,l,t)*phi(2,l,t)
      s2=(sw(2,1,l,t).hdot.phi(1,l,t))+sw(3,1,l,t)*phi(2,l,t)
      s3= sw(1,2,l,t)     *phi(3,l,t) +sw(2,2,l,t)*phi(4,l,t)
      s4=(sw(2,2,l,t).hdot.phi(3,l,t))+sw(3,2,l,t)*phi(4,l,t)
    */
    
    _vector_tensor_vector(v1,r->s0,s->s0);
    _vector_tensor_vector(v2,r->s0,s->s1);
    _vector_tensor_vector(v3,s->s0,r->s1);
    _su3_acc(v2,v3);
    _vector_tensor_vector(v3,r->s1,s->s1);
    
    _vector_tensor_vector(u1,r->s2,s->s2);
    _vector_tensor_vector(u2,r->s2,s->s3);
    _vector_tensor_vector(u3,s->s2,r->s3);
    _su3_acc(u2,u3);
    _vector_tensor_vector(u3,r->s3,s->s3);
    
    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0],u1,v1);
    _su3_plus_su3(lswp[1],u2,v2);
    _su3_plus_su3(lswp[2],u3,v3);
    _su3_dagger(lswp[3],lswp[1]);
    _su3_zero(lswp[1]);
    _su3_minus_su3(lswm[0],u1,v1);
    _su3_minus_su3(lswm[1],u2,v2);
    _su3_minus_su3(lswm[2],u3,v3);
    _su3_dagger(lswm[3],lswm[1]);
    _su3_zero(lswm[1]);
    
    /* add up the swm[0] and swp[0] */
    _su3_acc(swm[x][0],lswm[0]);
    _su3_acc(swm[x][1],lswm[1]);
    _su3_acc(swm[x][2],lswm[2]);
    _su3_acc(swm[x][3],lswm[3]);
    _su3_acc(swp[x][0],lswp[0]);
    _su3_acc(swp[x][1],lswp[1]);
    _su3_acc(swp[x][2],lswp[2]);
    _su3_acc(swp[x][3],lswp[3]);
  }
  return;
}

void sw_all(hamiltonian_field_t * const hf, const double kappa, const double c_sw) {
  int k,l;
  int x,xpk,xpl,xmk,xml,xpkml,xplmk,xmkml;
  su3 *w1,*w2,*w3,*w4;
  double ka_csw_8 = kappa*c_sw/8.;
  static su3 v1,v2,vv1,vv2,plaq;
  static su3 vis[4][4];
  
  for(x = 0; x < VOLUME; x++) {
    _minus_itimes_su3_plus_su3(vis[0][1],swm[x][1],swm[x][3]);
    _minus_su3_plus_su3(vis[0][2],swm[x][1],swm[x][3]);
    _itimes_su3_minus_su3(vis[0][3],swm[x][2],swm[x][0]);
    
    _minus_itimes_su3_plus_su3(vis[2][3],swp[x][1],swp[x][3]);
    _su3_plus_su3(vis[1][3],swp[x][1],swp[x][3]);
    
    _itimes_su3_minus_su3(vis[1][2],swp[x][2],swp[x][0]);

    /* anti-hermiticity */
    _su3_dagger(v1,vis[0][1]); _su3_minus_su3(vis[0][1],vis[0][1],v1);
    _su3_dagger(v1,vis[0][2]); _su3_minus_su3(vis[0][2],vis[0][2],v1);
    _su3_dagger(v1,vis[0][3]); _su3_minus_su3(vis[0][3],vis[0][3],v1);
    _su3_dagger(v1,vis[2][3]); _su3_minus_su3(vis[2][3],vis[2][3],v1);
    _su3_dagger(v1,vis[1][3]); _su3_minus_su3(vis[1][3],vis[1][3],v1);
    _su3_dagger(v1,vis[1][2]); _su3_minus_su3(vis[1][2],vis[1][2],v1);
    
    for(k = 0; k < 4; k++) {
      for(l=k+1;l<4;l++) {
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
 	_trace_lambda_mul_add_assign(hf->derivative[x][k], -ka_csw_8, vv1);

	_su3d_times_su3(vv2,*w1,vv1); 
	_su3_times_su3(vv1,vv2,*w1);
 	_trace_lambda_mul_add_assign(hf->derivative[xpk][l], -ka_csw_8, vv1);
	
	_su3_times_su3(vv2,vis[k][l],plaq); 
	_su3_dagger(vv1,vv2);
 	_trace_lambda_mul_add_assign(hf->derivative[x][l], -ka_csw_8, vv1);

	_su3d_times_su3(vv2,*w4,vv1); 
	_su3_times_su3(vv1,vv2,*w4);
 	_trace_lambda_mul_add_assign(hf->derivative[xpl][k], -ka_csw_8, vv1);
	
	w1=&hf->gaugefield[x][l];
	w2=&hf->gaugefield[xplmk][k];   /*dag*/
	w3=&hf->gaugefield[xmk][l];     /*dag*/
	w4=&hf->gaugefield[xmk][k];
	_su3_times_su3d(v1,*w1,*w2);
	_su3d_times_su3(v2,*w3,*w4);
	_su3_times_su3(plaq,v1,v2);
	
	_su3_times_su3(vv1,plaq,vis[k][l]);
 	_trace_lambda_mul_add_assign(hf->derivative[x][l], -ka_csw_8, vv1);
	
	_su3_dagger(vv1,v1); 
	_su3_times_su3d(vv2,vv1,vis[k][l]);
	_su3_times_su3d(vv1,vv2,v2);
 	_trace_lambda_mul_add_assign(hf->derivative[xplmk][k], -ka_csw_8, vv1);

	_su3_times_su3(vv2,*w3,vv1); 
	_su3_times_su3d(vv1,vv2,*w3);
 	_trace_lambda_mul_add_assign(hf->derivative[xmk][l], -ka_csw_8, vv1);

	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign(hf->derivative[xmk][k], -ka_csw_8, vv2);
	
	w1=&hf->gaugefield[xmk][k];   /*dag*/
	w2=&hf->gaugefield[xmkml][l]; /*dag*/
	w3=&hf->gaugefield[xmkml][k];
	w4=&hf->gaugefield[xml][l];
	_su3_times_su3(v1,*w2,*w1);
	_su3_times_su3(v2,*w3,*w4);
	
	_su3_times_su3d(vv1,*w1,vis[k][l]);
	_su3_times_su3d(vv2,vv1,v2);
	_su3_times_su3(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign(hf->derivative[xmk][k], -ka_csw_8, vv1);

	_su3_times_su3(vv2,*w2,vv1); 
	_su3_times_su3d(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign(hf->derivative[xmkml][l], -ka_csw_8, vv1);

	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign(hf->derivative[xmkml][k], -ka_csw_8, vv2);

	_su3d_times_su3(vv1,*w3,vv2); 
	_su3_times_su3(vv2,vv1,*w3);
 	_trace_lambda_mul_add_assign(hf->derivative[xml][l], -ka_csw_8, vv2);
	
	w1=&hf->gaugefield[xml][l];   /*dag*/
	w2=&hf->gaugefield[xml][k];
	w3=&hf->gaugefield[xpkml][l];
	w4=&hf->gaugefield[x][k];     /*dag*/
	_su3d_times_su3(v1,*w1,*w2);
	_su3_times_su3d(v2,*w3,*w4);
	
	_su3_times_su3d(vv1,*w1,vis[k][l]);
	_su3_times_su3d(vv2,vv1,v2);
	_su3_times_su3d(vv1,vv2,*w2);
 	_trace_lambda_mul_add_assign(hf->derivative[xml][l], -ka_csw_8, vv1);
	
	_su3_dagger(vv2,vv1);
 	_trace_lambda_mul_add_assign(hf->derivative[xml][k], -ka_csw_8, vv2);

	_su3d_times_su3(vv1,*w2,vv2); 
	_su3_times_su3(vv2,vv1,*w2);
 	_trace_lambda_mul_add_assign(hf->derivative[xpkml][l], -ka_csw_8, vv2);

	_su3_dagger(vv2,v2);  
	_su3_times_su3d(vv1,vv2,v1);
	_su3_times_su3d(vv2,vv1,vis[k][l]);
 	_trace_lambda_mul_add_assign(hf->derivative[x][k], -ka_csw_8, vv2);
      }
    }
  }
  return;
}

su3 * _swp;

int init_swpm(const int V) {
  int i=0;
  static int swpm_init=0;

  if(!swpm_init) {
    if((void*)(swp = (su3**)calloc(V, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(swm = (su3**)calloc(V, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(_swp = (su3*)calloc(2*4*V+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(2);
    }
#if (defined SSE || defined SSE2 || defined SSE3)
    swp[0] = (su3*)(((unsigned long int)(_swp)+ALIGN_BASE)&~ALIGN_BASE);
#else
    swp[0] = _swp;
#endif
    swm[0] = swp[0] + 4*V;
    for(i = 1; i < V; i++){
      swp[i] = swp[i-1]+4;
      swm[i] = swm[i-1]+4;
    }
    swpm_init = 1;
  }
  return(0);
}
