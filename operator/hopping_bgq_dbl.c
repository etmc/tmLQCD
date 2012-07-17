/**********************************************************************
 *
 *
 * Copyright (C) 2012 Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * and modified and extended by Carsten Urbach from 2003-2008
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
 *
 **********************************************************************/

#include "bgq.h"
#include "xlc_prefetch.h"

#define _hop_t_p()							\
  _vec_load2(r0, r1, r2, sp->s0);					\
  _vec_load2(r3, r4, r5, sp->s1);					\
  _vec_load2(r6, r7, r8, sp->s2);					\
  _vec_load2(r9, r10, r11, sp->s3);					\
  _vec_add_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11);	\
  _vec_su3_multiply_double2(up);					\
  _vec_cmplx_mul_double2(rs0, rs1, rs2, rs3, rs4, rs5, r6, r7, r8, r9, r10, r11, U0, ka0); \
  rs6 = rs0; rs7 = rs1; rs8 = rs2;					\
  rs9 = rs3; rs10= rs4; rs11= rs5;

#define _hop_t_m()							\
  _vec_load2(r0, r1, r2, sm->s0);					\
  _vec_load2(r3, r4, r5, sm->s1);					\
  _vec_load2(r6, r7, r8, sm->s2);					\
  _vec_load2(r9, r10, r11, sm->s3);					\
  _vec_sub_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11);	\
  _vec_su3_inverse_multiply_double2(um);				\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka0); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_sub_double2(rs6, rs7, rs8, rs9, rs10, rs11, r0, r1, r2, r3, r4, r5);

#define _hop_x_p()							\
  _vec_load2(r0, r1, r2, sp->s0);					\
  _vec_load2(r3, r4, r5, sp->s1);					\
  _vec_load2(r9, r10, r11, sp->s2);					\
  _vec_load2(r6, r7, r8, sp->s3);					\
  _vec_i_mul_add_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0); \
  _vec_su3_multiply_double2(up);					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka1); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r0, r1, r2, U1);

#define _hop_x_m()							\
  _vec_load2(r0, r1, r2, sm->s0);					\
  _vec_load2(r3, r4, r5, sm->s1);					\
  _vec_load2(r9, r10, r11, sm->s2);					\
  _vec_load2(r6, r7, r8, sm->s3);					\
  _vec_i_mul_sub_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0); \
  _vec_su3_inverse_multiply_double2(um);				\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka1); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r3, r4, r5, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r0, r1, r2, U1);

#define _hop_y_p()							\
  _vec_load2(r0, r1, r2, sp->s0);					\
  _vec_load2(r3, r4, r5, sp->s1);					\
  _vec_load2(r9, r10, r11, sp->s2);					\
  _vec_load2(r6, r7, r8, sp->s3);					\
  _vec_add2(r0, r1, r2, r6, r7, r8);					\
  _vec_sub2(r3, r4, r5, r9, r10, r11);					\
  _vec_su3_multiply_double2(up);					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka2); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5);	\
  _vec_sub2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_add2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_y_m()							\
  _vec_load2(r0, r1, r2, sm->s0);					\
  _vec_load2(r3, r4, r5, sm->s1);					\
  _vec_load2(r9, r10, r11, sm->s2);					\
  _vec_load2(r6, r7, r8, sm->s3);					\
  _vec_sub2(r0, r1, r2, r6, r7, r8);					\
  _vec_add2(r3, r4, r5, r9, r10, r11);					\
  _vec_su3_inverse_multiply_double2(um);				\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka2); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_add2(rs6, rs7, rs8, r3, r4, r5);					\
  _vec_sub2(rs9, rs10, rs11, r0, r1, r2);

#define _hop_z_p()							\
  _vec_load2(r0, r1, r2, sp->s0);					\
  _vec_load2(r3, r4, r5, sp->s1);					\
  _vec_load2(r6, r7, r8, sp->s2);					\
  _vec_load2(r9, r10, r11, sp->s3);					\
  _vec_i_mul_sub2(r3, r4, r5, r9, r10, r11, U1);			\
  _vec_i_mul_add2(r0, r1, r2, r6, r7, r8, U0);				\
  _vec_su3_multiply_double2(up);					\
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka3); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_sub2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_add2(rs9, rs10, rs11, r3, r4, r5, U1);

#define _hop_z_m()							\
  _vec_load2(r0, r1, r2, sm->s0);					\
  _vec_load2(r3, r4, r5, sm->s1);					\
  _vec_load2(r6, r7, r8, sm->s2);					\
  _vec_load2(r9, r10, r11, sm->s3);					\
  _vec_i_mul_sub2(r0, r1, r2, r6, r7, r8, U0);				\
  _vec_i_mul_add2(r3, r4, r5, r9, r10, r11, U1);			\
  _vec_su3_inverse_multiply_double2(um);				\
  _vec_cmplxcg_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, U0, ka3); \
  _vec_add_double2(rs0, rs1, rs2, rs3, rs4, rs5, r0, r1, r2, r3, r4, r5); \
  _vec_i_mul_add2(rs6, rs7, rs8, r0, r1, r2, U0);			\
  _vec_i_mul_sub2(rs9, rs10, rs11, r3, r4, r5, U1);


void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {

  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 * restrict ALIGN up;
  su3 * restrict ALIGN um;
  spinor * restrict ALIGN sp;
  spinor * restrict ALIGN sm;
  spinor * restrict ALIGN rn;
  /* We have 32 registers available */
  vector4double ALIGN r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;
  vector4double ALIGN rs0, rs1, rs2, rs3, rs4, rs5, rs6, rs7, rs8, rs9, rs10, rs11;
  vector4double ALIGN U0, U1, U2, U3, U4, U6, U7;
  vector4double rtmp;
#pragma disjoint(*sp, *sm, *rn, *up, *um, *l, *k)

  __alignx(16,l);
  __alignx(16,k);

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }
  ioff2 = (VOLUME+RAND)/2-ioff;

  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; 
  icy=g_lexic2eosub[iy];

  sp=k+icy;

#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif
  /**************** loop over all lattice sites ******************/
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){
    rn=l+(icx-ioff);
    ix=g_eo2lexic[icx];
    /*********************** direction +t ************************/
    iy=g_idn[ix][0]; 
    icy=g_lexic2eosub[iy];
#    if (!defined _GAUGE_COPY)
    um=&g_gauge_field[iy][0]; 
#    else
    um=up+1;
#    endif
    sm=k+icy;

    _hop_t_p();

    /*********************** direction -t ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;
    
    _hop_t_m();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#    else
    um = up+1;
#    endif
    sm=k+icy;

    _hop_x_p();

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;

    _hop_x_m();

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#    else
    um= up+1;
#    endif
    sm=k+icy;

    _hop_y_p();

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;

    _hop_y_m();

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#    else
    um=up+1;
#    endif
    sm=k+icy;

    _hop_z_p();

    /*********************** direction -3 ************************/

    icz=icx+1;
    if(icz==((VOLUME+RAND)/2+ioff)) icz=ioff;
    iz=g_eo2lexic[icz];
    iy=g_iup[iz][0]; icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up=&g_gauge_field[iz][0];
#    endif
    sp=k+icy;

    _hop_z_m();

    _vec_store2(rn->s0, rs0, rs1, rs2);
    _vec_store2(rn->s1, rs3, rs4, rs5);
    _vec_store2(rn->s2, rs6, rs7, rs8);
    _vec_store2(rn->s3, rs9, rs10, rs11);
  }
  return;
}
