/**********************************************************************
 *
 * $Id$
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


/* input on k; output on l */


void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  int ix, i;
  su3 * restrict U ALIGN;
  static spinor rs;
  spinor * restrict s ALIGN;
  halfspinor ** phi ALIGN;
#if defined OPTERON
  const int predist=2;
#else
  const int predist=1;
#endif
#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge();
  }
#endif
  /* We will run through the source vector now */
  /* instead of the solution vector            */
  s = k;
  _prefetch_spinor(s);

  if(ieo == 0) {
    U = g_gauge_field_copy[0][0];
  }
  else {
    U = g_gauge_field_copy[1][0];
  }
  phi = NBPointer[ieo];

  _prefetch_su3(U);
  /**************** loop over all lattice sites ******************/
  ix=0;
  for(i = 0; i < (VOLUME)/2; i++){

    /*********************** direction +0 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka0);
    _sse_store_nt_up((*phi[ix]).s1);
    U++;
    ix++;
    /*********************** direction -0 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +1 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    /*next not needed?*/
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka1);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -1 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +2 ************************/
    _prefetch_su3(U+predist);

    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka2);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;
    /*********************** direction -2 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s3);
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s2);
    _sse_vector_add();
    _sse_store_nt((*phi[ix]).s1);
    ix++;

    /*********************** direction +3 ************************/
    _prefetch_su3(U+predist);
    _prefetch_spinor(s+1);

    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*U));
    _sse_vector_cmplx_mul(ka3);
    _sse_store_nt_up((*phi[ix]).s1);
    ix++;
    U++;

    /*********************** direction -3 ************************/
    _sse_load((*s).s0);
    _sse_load_up((*s).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
    _sse_store_nt((*phi[ix]).s0);

    _sse_load((*s).s1);
    _sse_load_up((*s).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
    _sse_store_nt((*phi[ix]).s1);
    ix++;
    s++;
  }

#    if (defined MPI && !defined _NO_COMM)
  xchange_halffield(); 
#    endif
  s = l;
  phi = NBPointer[2 + ieo];
  if(ieo == 0) {
    U = g_gauge_field_copy[1][0];
  }
  else {
    U = g_gauge_field_copy[0][0];
  }
  _prefetch_su3(U);
  
  /* Now we sum up and expand to a full spinor */
  ix = 0;
  for(i = 0; i < (VOLUME)/2; i++){
    /*********************** direction +0 ************************/
    _vector_assign(rs.s0, (*phi[ix]).s0);
    _vector_assign(rs.s2, (*phi[ix]).s0);
    _vector_assign(rs.s1, (*phi[ix]).s1);
    _vector_assign(rs.s3, (*phi[ix]).s1);
    ix++;

    /*********************** direction -0 ************************/
    _prefetch_su3(U+predist);
      
    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);

    _sse_load((*phi[ix]).s1);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka0);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);

    ix++;
    U++;
    /*********************** direction +1 ************************/
    _sse_load_up((*phi[ix]).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s3); 

    _sse_load_up((*phi[ix]).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);       
    ix++;

    /*********************** direction -1 ************************/

    _prefetch_su3(U+predist);

    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka1);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka1);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s2);
    ix++;
    U++;

    /*********************** direction +2 ************************/
    _sse_load_up((*phi[ix]).s0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3);

    _sse_load_up((*phi[ix]).s1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);      
    ix++;

    /*********************** direction -2 ************************/

    _prefetch_su3(U+predist);

    _sse_load((*phi[ix]).s0);
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka2);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka2);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);      
    ix++;
    U++;
    /*********************** direction +3 ************************/
    _sse_load_up((*phi[ix]).s0);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);

    _sse_load_up((*phi[ix]).s1);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    ix++;
    /*********************** direction -3 ************************/

    _prefetch_su3(U+predist); 
    _prefetch_spinor(s+1);

    _sse_load((*phi[ix]).s0);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store_nt((*s).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store_nt((*s).s2);

    _sse_load((*phi[ix]).s1);
      
    _sse_su3_inverse_multiply((*U));
    _sse_vector_cmplxcg_mul(ka3);

    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store_nt((*s).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store_nt((*s).s3);
    ix++;
    U++;
    s++;
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}

