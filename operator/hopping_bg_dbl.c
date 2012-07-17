/**********************************************************************
 *
 *
 * Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 Carsten Urbach
 *
 * BG and halfspinor versions (C) 2007, 2008 Carsten Urbach
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


void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;
  spinor * restrict rn ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;

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
    /*********************** direction +0 ************************/
    iy=g_idn[ix][0]; 
    icy=g_lexic2eosub[iy];
#    if (!defined _GAUGE_COPY)
    um=&g_gauge_field[iy][0]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg1(sp->s1);
    _bgl_load_reg0_up(sp->s2);
    _bgl_load_reg1_up(sp->s3);
    _bgl_vector_add_reg0();
    _bgl_vector_add_reg1();
    /* result is now in regx0, regx1, regx2 x = 0,1 */

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka0);
    _bgl_store_reg0_up_rs0();
    _bgl_store_reg0_up_rs2();
    _bgl_store_reg1_up_rs1();
    _bgl_store_reg1_up_rs3();

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg1(sm->s1);
    _bgl_load_reg0_up(sm->s2);
    _bgl_load_reg1_up(sm->s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka0);

    _bgl_add_to_rs0_reg0();
    _bgl_sub_from_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs3_reg1();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#    else
    um = up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg1(sp->s1);
    _bgl_load_reg0_up(sp->s3);
    _bgl_load_reg1_up(sp->s2);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_add_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_i_mul_sub_from_rs3_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg1();

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg1(sm->s1);
    _bgl_load_reg0_up(sm->s3);
    _bgl_load_reg1_up(sm->s2);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_sub_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka1);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_add_to_rs3_reg0();
    _bgl_i_mul_add_to_rs2_reg1();      

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#    else
    um= up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg1(sp->s1);
    _bgl_load_reg1_up(sp->s2);
    _bgl_load_reg0_up(sp->s3);
    _bgl_vector_add_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs2_reg1();
    _bgl_add_to_rs3_reg0();


    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg1(sm->s1);
    _bgl_load_reg1_up(sm->s2);
    _bgl_load_reg0_up(sm->s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka2);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs2_reg1();
    _bgl_sub_from_rs3_reg0();

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#    else
    um=up+1;
#    endif
    _prefetch_su3(um); 
    sm=k+icy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0(sp->s0);
    _bgl_load_reg1(sp->s1);
    _bgl_load_reg0_up(sp->s2);
    _bgl_load_reg1_up(sp->s3);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka3);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg0();
    _bgl_i_mul_add_to_rs3_reg1();

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
    _prefetch_su3(up); 
    sp=k+icy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0(sm->s0);
    _bgl_load_reg1(sm->s1);
    _bgl_load_reg0_up(sm->s2);
    _bgl_load_reg1_up(sm->s3);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka3);


      
    _bgl_add_to_rs0_reg0();
    _bgl_store_rs0(rn->s0);
    _bgl_i_mul_add_to_rs2_reg0();
    _bgl_store_rs2(rn->s2);

    _bgl_add_to_rs1_reg1();
    _bgl_store_rs1(rn->s1);
    _bgl_i_mul_sub_from_rs3_reg1();
    _bgl_store_rs3(rn->s3);

    /************************ end of loop ************************/
  }
}
