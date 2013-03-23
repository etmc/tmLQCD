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


#    if (defined _USE_TSPLITPAR) /* needs also SSE */

/***********************************
 * 
 * Aurora version
 * Author: Luigi Scorzato (scorzato@ect.it)
 * (last modified 20.4.2009)
 * The strategy of the code is explained in the file Strategy.txt
 *
 ************************************/

/* 4. */
/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int icx,icz,ioff;
  int ix,iz;
  int x0,icx0,jj;
  su3 *restrict up;
  su3 * restrict um;
  spinor * restrict sp;
  spinor * restrict sm;
  spinor * restrict rn;

# if (defined MPI)
#  ifdef PARALLELX
#   define  REQC 4
#  elif defined PARALLELXY
#   define  REQC 8
#  elif defined PARALLELXYZ
#   define  REQC 12
#  endif
  MPI_Request requests[REQC];
  MPI_Status status[REQC];
# endif

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

  if(ieo == 0){ /* even out - odd in */
    ioff = 0;
  } 
  else{ /* odd out - even in */
    ioff = (VOLUME+RAND)/2;
  }

  /* Loop over time direction. This is the outmost loop */
  for(x0=0;x0<T;x0++){

    /* start the communication of the timslice borders (non-blocking send and receive)*/
#    if (defined MPI && !defined _NO_COMM)
   xchange_field_open(k, ieo, x0, requests, status);
#    endif
    

  /* loop over timeslice. At: contribution of timelike links  */
   icx0=g_1st_eot[x0][ieo];
   jj =0;
   um=&g_gauge_field_copyt[icx0][0]-1; /* allowed? */
   for(icx = icx0; icx < icx0+TEOSLICE; icx++){
     rn=l+(icx-ioff);
    /*********************** direction +0 ************************/

    sp=k+g_iup_eo[icx][0]; /* all sp,sm,up,um could be moved up */
    up=um+1;

    _sse_load(sp->s0);
    _sse_load_up(sp->s2);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rn->s0);
    _sse_store_up(rn->s2);      
      
    _sse_load(sp->s1);
    _sse_load_up(sp->s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_up(rn->s1);
    _sse_store_up(rn->s3); 

    /*********************** direction -0 ************************/

    sm=k+g_idn_eo[icx][0];
    um=up+1;

    _sse_load(sm->s0);
    _sse_load_up(sm->s2);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s2);
    _sse_vector_sub();
    _sse_store(rn->s2);
      
    _sse_load(sm->s1);
    _sse_load_up(sm->s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka0);
      
    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s3);
    _sse_vector_sub();
    _sse_store(rn->s3);
    jj++;
   } /* end of loop over timeslice (At)*/

       
   /* complete the communication of the timslice borders (and wait) */
#if (defined MPI && !defined _NO_COMM)
   xchange_field_close(requests, status, REQC); /*    MPI_Waitall */
#endif

   /* loop over timeslice. Bt: contribution of spacelike links  */
   um=&g_gauge_field_copys[icx0][0]-1;
   for(icx = icx0; icx < icx0+TEOSLICE; icx++){
    ix=g_eo2lexic[icx];
    rn=l+(icx-ioff);
    /*********************** direction +1 ************************/

    sp=k+g_iup_eo[icx][1];
    up=um+1;

    _sse_load(sp->s0);
    _sse_load_up(sp->s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rn->s3); 
      
    _sse_load(sp->s1);
    _sse_load_up(sp->s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka1);

    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rn->s2);       

    /*********************** direction -1 ************************/

    sm=k+g_idn_eo[icx][1];
    um=up+1;

    _sse_load(sm->s0);
    _sse_load_up(sm->s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rn->s3);

    _sse_load(sm->s1);
    _sse_load_up(sm->s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka1);
      
    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rn->s2);

    /*********************** direction +2 ************************/

    sp=k+g_iup_eo[icx][2];
    up=um+1;

    _sse_load(sp->s0);
    _sse_load_up(sp->s3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s3);
    _sse_vector_add();
    _sse_store(rn->s3);
      
    _sse_load(sp->s1);
    _sse_load_up(sp->s2);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka2);

    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s2);
    _sse_vector_sub();
    _sse_store(rn->s2);      

    /*********************** direction -2 ************************/

    sm=k+g_idn_eo[icx][2];
    um=up+1;

    _sse_load(sm->s0);
    _sse_load_up(sm->s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s3);
    _sse_vector_sub();
    _sse_store(rn->s3);
      
    _sse_load(sm->s1);
    _sse_load_up(sm->s2);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka2);
      
    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s2);
    _sse_vector_add();
    _sse_store(rn->s2);      
      
    /*********************** direction +3 ************************/

    sp=k+g_iup_eo[icx][3];
    up=um+1;

    _sse_load(sp->s0);
    _sse_load_up(sp->s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rn->s2);
      
    _sse_load(sp->s1);
    _sse_load_up(sp->s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(ka3);

    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rn->s3);
      
    /*********************** direction -3 ************************/

    sm=k+g_idn_eo[icx][3];
    um=up+1;

    _sse_load(sm->s0);
    _sse_load_up(sm->s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rn->s0);
    _sse_vector_add();
    _sse_store(rn->s0);

    _sse_load(rn->s2);
    _sse_vector_i_mul();
    _sse_vector_add();
    _sse_store(rn->s2);

    _sse_load(sm->s1);
    _sse_load_up(sm->s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rn->s1);
    _sse_vector_add();
    _sse_store(rn->s1);

    _sse_load(rn->s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rn->s3);
   }  /* end of loop over timeslice (Bt)*/
  } /* x0=0; x0<T */
}

#    endif
