#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "sw.h"

/* 
 * This function sets the array sw[][][]
 */

void sw_term(){

  int k,l;
  int x, xpk, xpl, xmk, xml, xpkml, xplmk, xmkml;
  su3 *w1, *w2, *w3, *w4;
  static su3 v1, v2, plaq;
  static su3 fkl[4][4];
  static su3 magnetic[4], electric[4];
  static su3 aux;

  /*  compute the clover-leave */
  /*  l  __   __
   *    |  | |  |
   *    |__| |__|
   *     __   __
   *    |  | |  |
   *    |__| |__| k  */

  for(x=0; x<VOLUME; x++){
    for(k=0;k<4;k++){
      for(l=k+1;l<4;l++){
	xpk=g_iup[x][k];
	xpl=g_iup[x][l];
	xmk=g_idn[x][k];
	xml=g_idn[x][l];
	xpkml=g_idn[xpk][l];
	xplmk=g_idn[xpl][k];
	xmkml=g_idn[xml][k];
	w1=&g_gauge_field[x][k];
	w2=&g_gauge_field[xpk][l];
	w3=&g_gauge_field[xpl][k];
	w4=&g_gauge_field[x][l];
	_su3_times_su3(v1, *w1, *w2);
	_su3_times_su3(v2, *w4, *w3);
	_su3_times_su3d(plaq, v1, v2);
	w1=&g_gauge_field[x][l];
	w2=&g_gauge_field[xplmk][k];
	w3=&g_gauge_field[xmk][l];
	w4=&g_gauge_field[xmk][k];
	_su3_times_su3d(v1, *w1, *w2);
	_su3d_times_su3(v2, *w3, *w4);
	_su3_times_su3_acc(plaq, v1, v2);
	w1=&g_gauge_field[xmk][k];
	w2=&g_gauge_field[xmkml][l];
	w3=&g_gauge_field[xmkml][k];
	w4=&g_gauge_field[xml][l];
	_su3_times_su3(v1, *w2, *w1);
	_su3_times_su3(v2, *w3, *w4);
	_su3d_times_su3_acc(plaq, v1, v2);
	w1=&g_gauge_field[xml][l];
	w2=&g_gauge_field[xml][k];
	w3=&g_gauge_field[xpkml][l];
	w4=&g_gauge_field[x][k];
	_su3d_times_su3(v1, *w1, *w2);
	_su3_times_su3d(v2, *w3, *w4);
	_su3_times_su3_acc(plaq, v1, v2);
	_su3_dagger(v2, plaq); 
	_su3_minus_su3(fkl[k][l], plaq, v2);
      }
    }

    _su3_one(sw[x][0][0]);
    _su3_one(sw[x][2][0]);
    _su3_one(sw[x][0][1]);
    _su3_one(sw[x][2][1]);

    for(k=1;k<4;k++){
      _su3_assign(electric[k], fkl[0][k]);
    }
    _su3_assign(magnetic[1], fkl[2][3]);
    _su3_minus_assign(magnetic[2], fkl[1][3]);
    _su3_assign(magnetic[3], fkl[1][2]);

    /*  upper left block matrix  */

    _itimes_su3_minus_su3(aux, electric[3], magnetic[3]);
    _su3_refac_acc(sw[x][0][0], g_ka_csw_8, aux);

    _itimes_su3_minus_su3(aux, electric[1], magnetic[1]);
    _su3_minus_su3(v2, electric[2], magnetic[2]); 
    _su3_acc(aux, v2);
    _real_times_su3(sw[x][1][0], g_ka_csw_8, aux);

    _itimes_su3_minus_su3(aux, magnetic[3], electric[3]);
    _su3_refac_acc(sw[x][2][0], g_ka_csw_8, aux);

    /*  lower right block matrix */

    _itimes_su3_plus_su3(aux, electric[3], magnetic[3]);
    _su3_refac_acc(sw[x][0][1], (-g_ka_csw_8), aux);

    _itimes_su3_plus_su3(aux, electric[1], magnetic[1]);
    _su3_plus_su3(v2, electric[2], magnetic[2]);
    _su3_acc(aux, v2);
    _real_times_su3(sw[x][1][1], (-g_ka_csw_8), aux);

    _itimes_su3_plus_su3(aux, magnetic[3], electric[3]);
    _su3_refac_acc(sw[x][2][1], g_ka_csw_8, aux);

  }
}

/*
  !--------------------------------------------------------------!
  !  The subroutine sw_invert is needed for the                  !
  !  even_odd preconditioned Dirac operator with SW improvement. !
  !  Details can be found in  the notes sw.ps on tsun.desy.de    !
  !  by P. Weisz and U. Wolff.                                   !
  !--------------------------------------------------------------!
  !  inversion in place of complex matrix a without pivoting     !
  !  triangularization by householder reflexions                 !
  !  inversion of triangular matrix                              !
  !  inverse reflexions                                          !
  !--------------------------------------------------------------!
  !  a square matrix,  dimensioned 0:n-1                          !
  !  itrouble is counted up,  when a dangerously small diagonal   !
  !  element is encountered in the tringular matrix              !
  !  has to be initialized outside                               !
  !                                                              !
  !  Author: U. Wolff,  adapted to fortran90 by S. Sint,  29/10/95 !
  !--------------------------------------------------------------!
  !  ported to C by M.Hasenbusch Wed Oct 24 15:46:46 MEST 2001   !
  !______________________________________________________________!
*/
#define nm1 5
int six_invert(complex a[6][6]){

  static complex d[nm1+1], u[nm1+1];
  static complex sigma, z;
  static double p[nm1+1];
  static double s, q;
  int i, j, k;
  int ifail;
  ifail=0;
  for(k=0;k<nm1;k++){
    s=0.0;
    for(j=k+1;j<=nm1;j++){
      s+=a[j][k].re*a[j][k].re+a[j][k].im*a[j][k].im; 
    }
    s=1.+s/(a[k][k].re*a[k][k].re+a[k][k].im*a[k][k].im);
    s=sqrt(s);
    sigma.re=s*a[k][k].re; 
    sigma.im=s*a[k][k].im;
    a[k][k].re+=sigma.re; 
    a[k][k].im+=sigma.im;
    p[k]=sigma.re*a[k][k].re+sigma.im*a[k][k].im;
    /* d[k]=-1./sigma */
    q=sigma.re*sigma.re+sigma.im*sigma.im;
    if(q<tiny_t){
      ifail++;
    }
    q=1./q;
    d[k].re=-q*sigma.re; 
    d[k].im= q*sigma.im; 
    /* reflect all columns to the right */
    for(j=k+1;j<=nm1;j++){
      z.re=0.; z.im=0.; 
      for(i=k;i<=nm1;i++){
	z.re+=a[i][k].re*a[i][j].re+a[i][k].im*a[i][j].im;
	z.im+=a[i][k].re*a[i][j].im-a[i][k].im*a[i][j].re;
      }
      q=1./p[k];
      z.re=q*z.re; 
      z.im=q*z.im; 
      for(i=k;i<=nm1;i++){
	a[i][j].re-=z.re*a[i][k].re-z.im*a[i][k].im;
	a[i][j].im-=z.re*a[i][k].im+z.im*a[i][k].re;
      }
    }
  }
  sigma.re=a[nm1][nm1].re;
  sigma.im=a[nm1][nm1].im;
  q=sigma.re*sigma.re+sigma.im*sigma.im;
  if(q<tiny_t){
    ifail++;
  }
  q=1./q;
  d[nm1].re= q*sigma.re; 
  d[nm1].im=-q*sigma.im; 
  /*  inversion of upper triangular matrix in place
      (diagonal elements done already): */

  for(k=nm1;k>=0;k--){
    for(i=k-1;i>=0;i--){
      z.re=0.; z.im=0.;
      for(j=i+1;j<k;j++)
	{
	  z.re+=a[i][j].re*a[j][k].re-a[i][j].im*a[j][k].im;
	  z.im+=a[i][j].re*a[j][k].im+a[i][j].im*a[j][k].re;
	}
      z.re+=a[i][k].re*d[k].re-a[i][k].im*d[k].im;
      z.im+=a[i][k].re*d[k].im+a[i][k].im*d[k].re;
      a[i][k].re=-z.re*d[i].re+z.im*d[i].im;
      a[i][k].im=-z.re*d[i].im-z.im*d[i].re;
    }
  }     
  /* execute reflexions in reverse order from the right: */

  a[nm1][nm1].re=d[nm1].re;
  a[nm1][nm1].im=d[nm1].im;
  for(k=nm1-1;k>=0;k--){
    for(j=k;j<=nm1;j++){
      u[j].re=a[j][k].re;
      u[j].im=a[j][k].im;
    }
    a[k][k].re=d[k].re;
    a[k][k].im=d[k].im;
    for(j=k+1;j<=nm1;j++){
      a[j][k].re=0.0;
      a[j][k].im=0.0;
    }
    for(i=0;i<=nm1;i++){
      z.re=0.0;
      z.im=0.0;
      for(j=k;j<=nm1;j++){
        z.re+=a[i][j].re*u[j].re-a[i][j].im*u[j].im;
        z.im+=a[i][j].re*u[j].im+a[i][j].im*u[j].re;
      }
      z.re=z.re/p[k];         /* normalization */
      z.im=z.im/p[k];         /* normalization */
      for(j=k;j<=nm1;j++){
        a[i][j].re-= z.re*u[j].re+z.im*u[j].im; /* reflexion */
        a[i][j].im-=-z.re*u[j].im+z.im*u[j].re; /* reflexion */
      }
    }
  }
  return ifail;
}
    
double six_det(complex a[6][6]){

  static complex sigma, z;
  static complex det;
  static double p[nm1+1];
  static double s, q, ds;
  int i, j, k;
  int ifail;
  ifail=0;
  /* compute the determinant:*/
  det.re=1.0;  det.im=0.0;

  for(k=0;k<nm1;k++){
    s=0.0;
    for(j=k+1;j<=nm1;j++){
      s+=a[j][k].re*a[j][k].re+a[j][k].im*a[j][k].im; 
    }
    s=1.+s/(a[k][k].re*a[k][k].re+a[k][k].im*a[k][k].im);
    s=sqrt(s);
    sigma.re=s*a[k][k].re; 
    sigma.im=s*a[k][k].im;

    /* determinant */
    ds = sigma.re*det.re-sigma.im*det.im;
    det.im = sigma.re*det.im+sigma.im*det.re;
    det.re = ds;
    q = sigma.re*sigma.re+sigma.im*sigma.im;
    if(q<tiny_t){
      ifail++;
    }
    a[k][k].re += sigma.re; 
    a[k][k].im += sigma.im;
    p[k] = sigma.re*a[k][k].re+sigma.im*a[k][k].im;

    /* reflect all columns to the right */
    for(j=k+1;j<=nm1;j++){
      z.re=0.; z.im=0.; 
      for(i = k; i <= nm1; i++){
	z.re+=a[i][k].re*a[i][j].re+a[i][k].im*a[i][j].im;
	z.im+=a[i][k].re*a[i][j].im-a[i][k].im*a[i][j].re;
      }
      q=1./p[k];
      z.re=q*z.re; 
      z.im=q*z.im; 
      for(i = k; i <= nm1; i++){
	a[i][j].re-=z.re*a[i][k].re-z.im*a[i][k].im;
	a[i][j].im-=z.re*a[i][k].im+z.im*a[i][k].re;
      }
    }
  }
  sigma.re=a[nm1][nm1].re;
  sigma.im=a[nm1][nm1].im;

  /* determinant */
  ds = sigma.re*det.re-sigma.im*det.im;
  det.im = sigma.re*det.im+sigma.im*det.re;
  det.re = ds;
  q = sigma.re*sigma.re+sigma.im*sigma.im;
  if(q < tiny_t){
    ifail++;
  }
  /*printf("%f %f \n", det.re, det.im);*/
  return det.re;
}

/*definitions needed for the functions sw_trace(int ieo) and sw_invert(int ieo)*/
#define _a_C(A, B, C) \
a[0+(A)][0+(B)]=(C).c00; \
a[0+(A)][1+(B)]=(C).c01; \
a[0+(A)][2+(B)]=(C).c02; \
a[1+(A)][0+(B)]=(C).c10; \
a[1+(A)][1+(B)]=(C).c11; \
a[1+(A)][2+(B)]=(C).c12; \
a[2+(A)][0+(B)]=(C).c20; \
a[2+(A)][1+(B)]=(C).c21; \
a[2+(A)][2+(B)]=(C).c22;

#define _C_a(A, B, C) \
(C).c00=a[0+(A)][0+(B)]; \
(C).c01=a[0+(A)][1+(B)]; \
(C).c02=a[0+(A)][2+(B)]; \
(C).c10=a[1+(A)][0+(B)]; \
(C).c11=a[1+(A)][1+(B)]; \
(C).c12=a[1+(A)][2+(B)]; \
(C).c20=a[2+(A)][0+(B)]; \
(C).c21=a[2+(A)][1+(B)]; \
(C).c22=a[2+(A)][2+(B)];

/*
 * Computes Tr log (1+T_{ee|oo})
 */
double sw_trace(int ieo){

  int i, x, icx, ioff;
  su3 *w;
  static su3 v2;
  static complex a[6][6];
  static double tra;
  static double ks, kc, tr, ts, tt;

  ks=0.0;
  kc=0.0;

  if(ieo==0){
    ioff=0;
  } 
  else
    {
      ioff=(VOLUME+RAND)/2;
    }
  for(icx=ioff;icx<(VOLUME/2+ioff);icx++){
    x=g_eo2lexic[icx];
    for(i=0;i<2;i++){
      w=&sw[x][0][i];     _a_C(0, 0, *w);
      w=&sw[x][1][i];     _a_C(0, 3, *w);
      _su3_dagger(v2, *w); _a_C(3, 0, v2);
      w=&sw[x][2][i];     _a_C(3, 3, *w);
      tra=log(six_det(a));

      tr=tra+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  }
  kc=ks+kc;
#ifdef MPI
  MPI_Allreduce(&kc,  &ks,  1,  MPI_DOUBLE,  MPI_SUM,  MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

/*
 * This computes the array
 * sw_inv[][][]
 */

void sw_invert(int ieo){

  int icx, ioff;
  int i, x;
  su3 *w;
  static su3 v2;
  static complex a[6][6];
  if(ieo==0){
    ioff=0;
  } 
  else{
    ioff=(VOLUME+RAND)/2;
  }
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    x=g_eo2lexic[icx];
    for(i=0;i<2;i++){
      w=&sw[x][0][i];     _a_C(0, 0, *w);
      w=&sw[x][1][i];     _a_C(0, 3, *w);
      _su3_dagger(v2, *w); _a_C(3, 0, v2);
      w=&sw[x][2][i];     _a_C(3, 3, *w);

      six_invert(a);

      /*  copy "a" back to sw_inv */
      w=&sw_inv[x][0][i]; _C_a(0, 0, *w);
      w=&sw_inv[x][1][i]; _C_a(0, 3, *w);
      w=&sw_inv[x][2][i]; _C_a(3, 3, *w);
    }
  }
}

/*
 * Here we compute the derivative of the
 * clover term
 */

void sw_deriv(int ieo){

  int ioff, icx;
  int x;
  static su3 lswp[4], lswm[4];

  /*  compute the clover-leave */
  /*  l  __   __
   *    |  | |  |
   *    |__| |__|
   *     __   __
   *    |  | |  |
   *    |__| |__| k  */


  /* convention: Tr colver-leaf times insertion */
  /* task : put the matrix in question to the front */
  if(ieo==0){
    ioff=0;
  } 
  else{
    ioff=(VOLUME+RAND)/2;
  }
  for(icx=ioff;icx<(VOLUME/2+ioff);icx++){
    x=g_eo2lexic[icx];
    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0], sw_inv[x][0][1], sw_inv[x][0][0]);
    _su3_plus_su3(lswp[1], sw_inv[x][1][1], sw_inv[x][1][0]);
    _su3_plus_su3(lswp[2], sw_inv[x][2][1], sw_inv[x][2][0]);
    _su3_dagger(lswp[3], lswp[1]);
    _su3_assign(lswp[1], lswp[3]);
    _su3_minus_su3(lswm[0], sw_inv[x][0][1], sw_inv[x][0][0]);
    _su3_minus_su3(lswm[1], sw_inv[x][1][1], sw_inv[x][1][0]);
    _su3_minus_su3(lswm[2], sw_inv[x][2][1], sw_inv[x][2][0]);
    _su3_dagger(lswm[3], lswm[1]);
    _su3_assign(lswm[1], lswm[3]);

    /* add up the swm[0]  and swp[2] */
    _su3_acc(swm[x][0], lswm[0]);
    _su3_acc(swm[x][1], lswm[1]);
    _su3_acc(swm[x][2], lswm[2]);
    _su3_acc(swm[x][3], lswm[3]);
    _su3_acc(swp[x][0], lswp[0]);
    _su3_acc(swp[x][1], lswp[1]);
    _su3_acc(swp[x][2], lswp[2]);
    _su3_acc(swp[x][3], lswp[3]);

  }
}


/* ieo : even =0 ; odd =1  k: spinor-field on the right; l: on the left */
void sw_spinor(int ieo,  int kk,  int ll){

  int ioff;
  int icx;
  int x;
  spinor *r, *s;
  static su3 v1, v2, v3;
  static su3 u1, u2, u3;
  static su3 lswp[4], lswm[4];

  /*  compute the clover-leave */
  /*  l  __   __
   *    |  | |  |
   *    |__| |__|
   *     __   __
   *    |  | |  |
   *    |__| |__| k  */


  /* convention: Tr colver-leaf times insertion */
  /* task : put the matrix in question to the front */

  if(ieo==0){
    ioff=0;
  } 
  else{
    ioff=(VOLUME+RAND)/2;
  }
  /**************** loop over half of the lattice sites **/


  for(icx=ioff;icx<(VOLUME/2+ioff);icx++){
    x=g_eo2lexic[icx];
    r=&spinor_field[kk][icx-ioff];
    s=&spinor_field[ll][icx-ioff];

    /*
      s1= sw(1, 1, l, t)     *phi(1, l, t) +sw(2, 1, l, t)*phi(2, l, t)
      s2=(sw(2, 1, l, t).hdot.phi(1, l, t))+sw(3, 1, l, t)*phi(2, l, t)
      s3= sw(1, 2, l, t)     *phi(3, l, t) +sw(2, 2, l, t)*phi(4, l, t)
      s4=(sw(2, 2, l, t).hdot.phi(3, l, t))+sw(3, 2, l, t)*phi(4, l, t)
    */

    _vector_tensor_vector(v1, (*r).s0, (*s).s0);
    _vector_tensor_vector(v2, (*r).s0, (*s).s1);
    _vector_tensor_vector(v3, (*s).s0, (*r).s1);
    _su3_acc(v2, v3);
    _vector_tensor_vector(v3, (*r).s1, (*s).s1);

    _vector_tensor_vector(u1, (*r).s2, (*s).s2);
    _vector_tensor_vector(u2, (*r).s2, (*s).s3);
    _vector_tensor_vector(u3, (*s).s2, (*r).s3);
    _su3_acc(u2, u3);
    _vector_tensor_vector(u3, (*r).s3, (*s).s3);

    /* compute the insertion matrix */
    _su3_plus_su3(lswp[0], u1, v1);
    _su3_plus_su3(lswp[1], u2, v2);
    _su3_plus_su3(lswp[2], u3, v3);
    _su3_dagger(lswp[3], lswp[1]);
    _su3_zero(lswp[1]);
    _su3_minus_su3(lswm[0], u1, v1);
    _su3_minus_su3(lswm[1], u2, v2);
    _su3_minus_su3(lswm[2], u3, v3);
    _su3_dagger(lswm[3], lswm[1]);
    _su3_zero(lswm[1]);

    /* add up the swm[0] and swp[0] */
    _su3_acc(swm[x][0], lswm[0]);
    _su3_acc(swm[x][1], lswm[1]);
    _su3_acc(swm[x][2], lswm[2]);
    _su3_acc(swm[x][3], lswm[3]);
    _su3_acc(swp[x][0], lswp[0]);
    _su3_acc(swp[x][1], lswp[1]);
    _su3_acc(swp[x][2], lswp[2]);
    _su3_acc(swp[x][3], lswp[3]);
  }
}

void sw_all(){

  int k, l;
  int x, xpk, xpl, xmk, xml, xpkml, xplmk, xmkml;
  su3 *w1, *w2, *w3, *w4;
  static su3 v1, v2, vv1, vv2, plaq;
  static su3 vis[4][4];
  static su3adj resu;
  su3adj *der;

  for(x=0;x<VOLUME;x++){
    _minus_itimes_su3_plus_su3(vis[0][1], swm[x][1], swm[x][3]);
    _minus_su3_plus_su3(vis[0][2], swm[x][1], swm[x][3]);

    _itimes_su3_minus_su3(vis[0][3], swm[x][2], swm[x][0]);

    _minus_itimes_su3_plus_su3(vis[2][3], swp[x][1], swp[x][3]);
    _su3_plus_su3(vis[1][3], swp[x][1], swp[x][3]);

    _itimes_su3_minus_su3(vis[1][2], swp[x][2], swp[x][0]);

    /* anti-hermiticity */
    _su3_dagger(v1, vis[0][1]); 
    _su3_minus_su3(vis[0][1], vis[0][1], v1);
    _su3_dagger(v1, vis[0][2]); 
    _su3_minus_su3(vis[0][2], vis[0][2], v1);
    _su3_dagger(v1, vis[0][3]); 
    _su3_minus_su3(vis[0][3], vis[0][3], v1);
    _su3_dagger(v1, vis[2][3]); 
    _su3_minus_su3(vis[2][3], vis[2][3], v1);
    _su3_dagger(v1, vis[1][3]); 
    _su3_minus_su3(vis[1][3], vis[1][3], v1);
    _su3_dagger(v1, vis[1][2]); 
    _su3_minus_su3(vis[1][2], vis[1][2], v1);

    for(k=0;k<4;k++){
      for(l=k+1;l<4;l++){
	xpk=g_iup[x][k];
	xpl=g_iup[x][l];
	xmk=g_idn[x][k];
	xml=g_idn[x][l];
	xpkml=g_idn[xpk][l];
	xplmk=g_idn[xpl][k];
	xmkml=g_idn[xml][k];
	w1=&g_gauge_field[x][k];
	w2=&g_gauge_field[xpk][l];
	w3=&g_gauge_field[xpl][k];   /*dag*/
	w4=&g_gauge_field[x][l];     /*dag*/

	_su3_times_su3(v1, *w1, *w2);
	_su3_times_su3(v2, *w4, *w3);
	_su3_times_su3d(plaq, v1, v2);

	_su3_times_su3(vv1, plaq, vis[k][l]);
	_trace_lambda(resu, vv1); 
	der=&dclover[x][k]; 
	_add_su3adj(*der, resu);

	_su3d_times_su3(vv2, *w1, vv1); 
	_su3_times_su3(vv1, vv2, *w1);
	_trace_lambda(resu, vv1); 
	der=&dclover[xpk][l]; 
	_add_su3adj(*der, resu);
 
	_su3_times_su3(vv2, vis[k][l], plaq); 
	_su3_dagger(vv1, vv2);
	_trace_lambda(resu, vv1); 
	der=&dclover[x][l]; 
	_add_su3adj(*der, resu);
	_su3d_times_su3(vv2, *w4, vv1); 
	_su3_times_su3(vv1, vv2, *w4);
	_trace_lambda(resu, vv1); 
	der=&dclover[xpl][k]; 
	_add_su3adj(*der, resu);

	w1=&g_gauge_field[x][l];
	w2=&g_gauge_field[xplmk][k];   /*dag*/
	w3=&g_gauge_field[xmk][l];     /*dag*/
	w4=&g_gauge_field[xmk][k];
	_su3_times_su3d(v1, *w1, *w2);
	_su3d_times_su3(v2, *w3, *w4);
	_su3_times_su3(plaq, v1, v2);

	_su3_times_su3(vv1, plaq, vis[k][l]);
	_trace_lambda(resu, vv1); 
	der=&dclover[x][l]; 
	_add_su3adj(*der, resu);

	_su3_dagger(vv1, v1); 
	_su3_times_su3d(vv2, vv1, vis[k][l]);
	_su3_times_su3d(vv1, vv2, v2);
	_trace_lambda(resu, vv1); 
	der=&dclover[xplmk][k]; 
	_add_su3adj(*der, resu);

	_su3_times_su3(vv2, *w3, vv1); 
	_su3_times_su3d(vv1, vv2, *w3);
	_trace_lambda(resu, vv1); 
	der=&dclover[xmk][l]; 
	_add_su3adj(*der, resu);

	_su3_dagger(vv2, vv1);
	_trace_lambda(resu, vv2); 
	der=&dclover[xmk][k]; 
	_add_su3adj(*der, resu);

	w1=&g_gauge_field[xmk][k];   /*dag*/
	w2=&g_gauge_field[xmkml][l]; /*dag*/
	w3=&g_gauge_field[xmkml][k];
	w4=&g_gauge_field[xml][l];
	_su3_times_su3(v1, *w2, *w1);
	_su3_times_su3(v2, *w3, *w4);

	_su3_times_su3d(vv1, *w1, vis[k][l]);
	_su3_times_su3d(vv2, vv1, v2);
	_su3_times_su3(vv1, vv2, *w2);
	_trace_lambda(resu, vv1); 
	der=&dclover[xmk][k]; 
	_add_su3adj(*der, resu);

	_su3_times_su3(vv2, *w2, vv1); 
	_su3_times_su3d(vv1, vv2, *w2);
	_trace_lambda(resu, vv1); 
	der=&dclover[xmkml][l]; 
	_add_su3adj(*der, resu);

	_su3_dagger(vv2, vv1);
	_trace_lambda(resu, vv2); 
	der=&dclover[xmkml][k]; 
	_add_su3adj(*der, resu);

	_su3d_times_su3(vv1, *w3, vv2); 
	_su3_times_su3(vv2, vv1, *w3);
	_trace_lambda(resu, vv2); 
	der=&dclover[xml][l]; 
	_add_su3adj(*der, resu);

	w1=&g_gauge_field[xml][l];   /*dag*/
	w2=&g_gauge_field[xml][k];
	w3=&g_gauge_field[xpkml][l];
	w4=&g_gauge_field[x][k];     /*dag*/
	_su3d_times_su3(v1, *w1, *w2);
	_su3_times_su3d(v2, *w3, *w4);

	_su3_times_su3d(vv1, *w1, vis[k][l]);
	_su3_times_su3d(vv2, vv1, v2);
	_su3_times_su3d(vv1, vv2, *w2);
	_trace_lambda(resu, vv1); 
	der=&dclover[xml][l]; 
	_add_su3adj(*der, resu);

	_su3_dagger(vv2, vv1);
	_trace_lambda(resu, vv2); 
	der=&dclover[xml][k]; 
	_add_su3adj(*der, resu);

	_su3d_times_su3(vv1, *w2, vv2); 
	_su3_times_su3(vv2, vv1, *w2);
	_trace_lambda(resu, vv2); 
	der=&dclover[xpkml][l]; 
	_add_su3adj(*der, resu);

	_su3_dagger(vv2, v2);  _su3_times_su3d(vv1, vv2, v1);
	_su3_times_su3d(vv2, vv1, vis[k][l]);
	_trace_lambda(resu, vv2); 
	der=&dclover[x][k]; 
	_add_su3adj(*der, resu);
      }
    }
  }
}
