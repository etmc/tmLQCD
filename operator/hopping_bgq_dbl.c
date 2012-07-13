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


{
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;
  spinor * restrict rn ALIGN;
  /* We have 32 registers available */
  vector4double r[12];
  vector4double U[9];
  /* The following contains the result spinor */
  vector4double rs[12];

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

    vec_load2(r, &sp->s0);
    vec_load2(r+3, &sp->s1);
    vec_load2(r+6, &sp->s2);
    vec_load2(r+9, &sp->s3);
    // s0 + s2 and s1 + s3
    vec_add_double2(r, &r[6]);
    // result is now in r[0-5] 
    vec_su3_multiply_double2(up, U, r);
    // result is now in r[6-11]
    // mult with ka0 and store in rs
    vec_cmplx_mul_double2(rs, &r[6], U, &ka0);
    rs[6] = rs[0]; rs[7] = rs[1]; rs[8] = rs[2];
    rs[9] = rs[3]; rs[10]= rs[4]; rs[11]= rs[5];

    /*********************** direction -t ************************/

    iy=g_iup[ix][1]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;
    vec_load2(r, &sm->s0);
    vec_load2(r+3, &sm->s1);
    vec_load2(r+6, &sm->s2);
    vec_load2(r+9, &sm->s3);
    // s0 - s2 and s1 - s3
    vec_sub_double2(r, &r[6]);
    // result is now in r[0-5]
    vec_su3_inverse_multiply_double2(um, U, r);
    // result is now in r[6-11]
    // mult with k0
    vec_cmplxcg_mul_double2(r, &r[6], U, &ka0);
    // result in r[0-5] now
    vec_add_double2(rs, r);
    vec_sub_double2(&rs[6], r);

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1]; 
#    else
    um = up+1;
#    endif
    sm=k+icy;

    vec_load2(r, &sp->s0);
    vec_load2(r+3, &sp->s1);
    vec_load2(r+9, &sp->s2);
    vec_load2(r+6, &sp->s3);
    vec_i_mul_add_double2(r, &r[6], U);
    vec_su3_multiply_double2(up, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka1);
    vec_add_double2(rs, r);
    vec_i_mul_sub2(&rs[6], &r[3], U);
    vec_i_mul_sub2(&rs[9], &r[0], U);

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;

    vec_load2(r, &sm->s0);
    vec_load2(r+3, &sm->s1);
    vec_load2(r+9, &sm->s2);
    vec_load2(r+6, &sm->s3);
    vec_i_mul_sub_double2(r, &r[6], U);
    vec_su3_inverse_multiply_double2(um, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka1);
    vec_add_double2(rs, r);
    vec_i_mul_add2(&rs[6], &r[3], U);
    vec_i_mul_add2(&rs[9], &r[0], U);

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][2]; 
#    else
    um= up+1;
#    endif
    sm=k+icy;

    vec_load2(r, &sp->s0);
    vec_load2(r+3, &sp->s1);
    vec_load2(r+9, &sp->s2);
    vec_load2(r+6, &sp->s3);
    vec_add2(r, &r[6]);
    vec_sub2(r+3, &r[9]);
    vec_su3_multiply_double2(up, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka2);
    vec_add_double2(rs, r);
    vec_sub2(&rs[6], &r[3]);
    vec_add2(&rs[9], r);

    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 
    icy=g_lexic2eosub[iy];

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif
    sp=k+icy;

    vec_load2(r, &sm->s0);
    vec_load2(r+3, &sm->s1);
    vec_load2(r+9, &sm->s2);
    vec_load2(r+6, &sm->s3);
    vec_sub2(r, r+6);
    vec_add2(r+3, r+9);
    vec_su3_inverse_multiply_double2(um, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka2);
    vec_add_double2(rs, r);
    vec_add2(rs+6, r+3);
    vec_sub2(rs+9, r);

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 
    icy=g_lexic2eosub[iy];

#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][3]; 
#    else
    um=up+1;
#    endif
    sm=k+icy;

    vec_load2(r, &sp->s0);
    vec_load2(r+3, &sp->s1);
    vec_load2(r+6, &sp->s2);
    vec_load2(r+9, &sp->s3);
    vec_i_mul_add2(r, r+6, U);
    vec_i_mul_sub2(r+3, r+9, U);
    vec_su3_multiply_double2(up, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka3);
    vec_add_double2(rs, r);
    vec_i_mul_sub2(rs+6, r, U);
    vec_i_mul_add2(rs+9, r+3, U);

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

    vec_load2(r, &sm->s0);
    vec_load2(r+3, &sm->s1);
    vec_load2(r+6, &sm->s2);
    vec_load2(r+9, &sm->s3);
    vec_i_mul_sub2(r, r+6, U);
    vec_i_mul_add2(r+3, r+9, U);
    vec_su3_inverse_multiply_double2(um, U, r);
    vec_cmplx_mul_double2(r, &r[6], U, &ka3);
    vec_add_double2(rs, r);
    vec_store2(&rn->s0, rs);
    vec_store2(&rn->s1, rs+3);
    vec_i_mul_add2(rs+6, r, U);
    vec_store2(&rn->s2, rs+6);
    vec_i_mul_sub2(rs+9, r+3, U);
    vec_store2(&rn->s3, rs+9);
  }
}
