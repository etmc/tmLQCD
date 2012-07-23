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


/* l output , k input*/
/* for ieo=0, k resides on  odd sites and l on even sites */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
  if(k == l){
    printf("Error in H_psi (simple.c):\n");
    printf("Arguments k and l must be different\n");
    printf("Program aborted\n");
    exit(1);
  }

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

#ifdef OMP
#define static
#pragma omp parallel
  {
  su3 * restrict U0 ALIGN;
#endif

  int i,ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  spinor rs;
  static su3_vector psi, chi, psi2, chi2;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;

#ifdef OMP
#undef static
#endif

#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif
#ifdef XLC
#pragma disjoint(*l, *k)
#pragma disjoint(*k, *U)
#pragma disjoint(*l, *U)
#pragma disjoint(*U, *s)
#pragma disjoint(*k, *s)
#pragma disjoint(*l, *s)
  __alignx(32, l);
  __alignx(32, k);
  __alignx(32, U);
  __alignx(32, s);
#endif

#ifndef OMP  
  s = k;
  if(ieo == 0) {
    U = g_gauge_field_copy[0][0];
  }
  else {
    U = g_gauge_field_copy[1][0];
  }
#else
  if(ieo == 0) {
    U0 = g_gauge_field_copy[0][0];
  }
  else {
    U0 = g_gauge_field_copy[1][0];
  }
#endif

  if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
    phi32 = NBPointer32[ieo];
      
    /**************** loop over all lattice sites ****************/
#ifdef OMP
#pragma omp for
#else
    ix=0;
#endif
    for(i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
      U=U0+i*4;
      s=k+i;
      ix=i*8;
#endif
      _vector_assign(rs.s0, s->s0);
      _vector_assign(rs.s1, s->s1);
      _vector_assign(rs.s2, s->s2);
      _vector_assign(rs.s3, s->s3);
      s++;
      /*********************** direction +0 ************************/
      
      _vector_add(psi, rs.s0, rs.s2);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector(phi32[ix]->s0, ka0, chi);
      
      _vector_add(psi, rs.s1, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector(phi32[ix]->s1, ka0, chi);
            
      U++;
      ix++;
    
      /*********************** direction -0 ************************/

      _vector_sub(phi32[ix]->s0, rs.s0, rs.s2);
      _vector_sub(phi32[ix]->s1, rs.s1, rs.s3);

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);

      _su3_multiply(chi, (*U), psi);
      _complex_times_vector(phi32[ix]->s0, ka1, chi);

      _vector_i_add(psi, rs.s1, rs.s2);

      _su3_multiply(chi, (*U), psi);
      _complex_times_vector(phi32[ix]->s1, ka1, chi);

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub(phi32[ix]->s0, rs.s0, rs.s3);
      _vector_i_sub(phi32[ix]->s1, rs.s1, rs.s2);

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector(phi32[ix]->s0, ka2, chi);

      _vector_sub(psi, rs.s1, rs.s2);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector(phi32[ix]->s1, ka2, chi);
      
      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub(phi32[ix]->s0, rs.s0, rs.s3);
      _vector_add(phi32[ix]->s1, rs.s1, rs.s2);
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);
      
      _su3_multiply(chi, (*U), psi);
      _complex_times_vector(phi32[ix]->s0, ka3, chi);


      _vector_i_sub(psi, rs.s1, rs.s3);

      _su3_multiply(chi,(*U),psi);
      _complex_times_vector(phi32[ix]->s1, ka3, chi);

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub(phi32[ix]->s0, rs.s0, rs.s2);
      _vector_i_add(phi32[ix]->s1, rs.s1, rs.s3);

#ifndef OMP
      ix++;
#endif
      /************************ end of loop ************************/
    }

#ifdef OMP
#pragma omp single
{
#endif

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield32(); 
#    endif

#ifdef OMP
}
#endif

#ifndef OMP
    s = l;
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
#else
    if(ieo == 0) {
      U0 = g_gauge_field_copy[1][0];
    }
    else {
      U0 = g_gauge_field_copy[0][0];
    }
#endif
   
   phi32 = NBPointer32[2 + ieo];

#ifdef OMP
#pragma omp for
#else
    ix = 0;
#endif
    for(i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
      ix=i*8;
      s=l+i;
      U=U0+i*4;
#endif
      /*********************** direction +0 ************************/
      _vector_assign(rs.s0, phi32[ix]->s0);
      _vector_assign(rs.s2, phi32[ix]->s0);
      _vector_assign(rs.s1, phi32[ix]->s1);
      _vector_assign(rs.s3, phi32[ix]->s1);
      ix++;
      /*********************** direction -0 ************************/
      _vector_assign(psi, phi32[ix]->s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka0,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s2, psi);

      _vector_assign(psi, phi32[ix]->s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka0,chi);
      
      _vector_add_assign(rs.s1, psi);
      _vector_sub_assign(rs.s3, psi);
      ix++;
      U++;
      /*********************** direction +1 ************************/

      _vector_add_assign(rs.s0, phi32[ix]->s0);
      _vector_i_sub_assign(rs.s3, phi32[ix]->s0);

      _vector_add_assign(rs.s1, phi32[ix]->s1);
      _vector_i_sub_assign(rs.s2, phi32[ix]->s1);
    
      ix++;
      /*********************** direction -1 ************************/
      _vector_assign(psi, phi32[ix]->s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka1,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_i_add_assign(rs.s3, psi);

      _vector_assign(psi, phi32[ix]->s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka1,chi);

      _vector_add_assign(rs.s1, psi);
      _vector_i_add_assign(rs.s2, psi);

      U++;
      ix++;

      /*********************** direction +2 ************************/

      _vector_add_assign(rs.s0, phi32[ix]->s0);
      _vector_add_assign(rs.s3, phi32[ix]->s0);

      _vector_add_assign(rs.s1, phi32[ix]->s1);
      _vector_sub_assign(rs.s2, phi32[ix]->s1);
    
      ix++;
      /*********************** direction -2 ************************/

      _vector_assign(psi, phi32[ix]->s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka2,chi);

      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s3, psi);

      _vector_assign(psi, phi32[ix]->s1);
      _su3_inverse_multiply(chi, (*U), psi);
      _complexcjg_times_vector(psi,ka2,chi);
      
      _vector_add_assign(rs.s1, psi);
      _vector_add_assign(rs.s2, psi);

      U++;
      ix++;
      /*********************** direction +3 ************************/

      _vector_add_assign(rs.s0, phi32[ix]->s0);
      _vector_i_sub_assign(rs.s2, phi32[ix]->s0);

      _vector_add_assign(rs.s1, phi32[ix]->s1);
      _vector_i_add_assign(rs.s3, phi32[ix]->s1);

      ix++;

      /*********************** direction -3 ************************/

      _vector_assign(psi, phi32[ix]->s0);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka3,chi);
      
      _vector_add(s->s0, rs.s0, psi);
      _vector_i_add(s->s2, rs.s2, psi);

      _vector_assign(psi, phi32[ix]->s1);
      _su3_inverse_multiply(chi,(*U), psi);
      _complexcjg_times_vector(psi,ka3,chi);

      _vector_add(s->s1, rs.s1, psi);
      _vector_i_sub(s->s3, rs.s3, psi);

#ifndef OMP
      U++;
      ix++;
      s++;
#endif
    }
  }
  else {
    phi = NBPointer[ieo];
      
    /**************** loop over all lattice sites ****************/
#ifdef OMP
#pragma omp for
#else
    ix=0;
#endif
    /* #pragma ivdep*/
    for(i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
      s=k+i;
      ix=i*8;
      U=U0+i*4;
#endif
      _vector_assign(rs.s0, s->s0);
      _vector_assign(rs.s1, s->s1);
      _vector_assign(rs.s2, s->s2);
      _vector_assign(rs.s3, s->s3);
      s++;
      /*********************** direction +0 ************************/
      
      _vector_add(psi, rs.s0, rs.s2);
      _vector_add(psi2, rs.s1, rs.s3);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector(phi[ix]->s0, ka0, chi);
      _complex_times_vector(phi[ix]->s1, ka0, chi2);
            
      U++;
      ix++;
    
      /*********************** direction -0 ************************/

      _vector_sub(phi[ix]->s0, rs.s0, rs.s2);
      _vector_sub(phi[ix]->s1, rs.s1, rs.s3);

      ix++;

      /*********************** direction +1 ************************/

      _vector_i_add(psi, rs.s0, rs.s3);
      _vector_i_add(psi2, rs.s1, rs.s2);
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2, (*U), psi2);
      _complex_times_vector(phi[ix]->s0, ka1, chi);
      _complex_times_vector(phi[ix]->s1, ka1, chi2);

      U++;
      ix++;

      /*********************** direction -1 ************************/

      _vector_i_sub(phi[ix]->s0, rs.s0, rs.s3);
      _vector_i_sub(phi[ix]->s1, rs.s1, rs.s2);

      ix++;
      /*********************** direction +2 ************************/

      _vector_add(psi, rs.s0, rs.s3);
      _vector_sub(psi2, rs.s1, rs.s2);
      _su3_multiply(chi,(*U),psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector(phi[ix]->s0, ka2, chi);
      _complex_times_vector(phi[ix]->s1, ka2, chi2);
      
      U++;
      ix++;

      /*********************** direction -2 ************************/

      _vector_sub(phi[ix]->s0, rs.s0, rs.s3);
      _vector_add(phi[ix]->s1, rs.s1, rs.s2);
      ix++;

      /*********************** direction +3 ************************/

      _vector_i_add(psi, rs.s0, rs.s2);
      _vector_i_sub(psi2, rs.s1, rs.s3);      
      _su3_multiply(chi, (*U), psi);
      _su3_multiply(chi2,(*U),psi2);
      _complex_times_vector(phi[ix]->s0, ka3, chi);
      _complex_times_vector(phi[ix]->s1, ka3, chi2);

      U++;
      ix++;
      /*********************** direction -3 ************************/

      _vector_i_sub(phi[ix]->s0, rs.s0, rs.s2);
      _vector_i_add(phi[ix]->s1, rs.s1, rs.s3);

#ifndef OMP
      ix++;
#endif
      /************************ end of loop ************************/
    }

#ifdef OMP
#pragma omp single
{
#endif

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield(); 
#    endif

#ifdef OMP
}
#endif

#ifndef OMP
    s = l;
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
#else
    if(ieo == 0) {
      U0 = g_gauge_field_copy[1][0];
    }
    else {
      U0 = g_gauge_field_copy[0][0];
    }
#endif

    phi = NBPointer[2 + ieo];

#ifdef OMP
#pragma omp for
#else
    ix = 0;
#endif
    /* #pragma ivdep */
    for(i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
      ix=i*8;
      U=U0+i*4;
      s=l+i;
#endif

      /*********************** direction +0 ************************/
      _vector_assign(rs.s0, phi[ix]->s0);
      _vector_assign(rs.s2, phi[ix]->s0);
      _vector_assign(rs.s1, phi[ix]->s1);
      _vector_assign(rs.s3, phi[ix]->s1);
      ix++;
      /*********************** direction -0 ************************/
      _su3_inverse_multiply(chi,(*U),phi[ix]->s0);
      _su3_inverse_multiply(chi2,(*U),phi[ix]->s1);
      _complexcjg_times_vector(psi,ka0,chi);
      _complexcjg_times_vector(psi2,ka0,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s2, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_sub_assign(rs.s3, psi2);
      ix++;
      U++;
      /*********************** direction +1 ************************/

      _vector_add_assign(rs.s0, phi[ix]->s0);
      _vector_i_sub_assign(rs.s3, phi[ix]->s0);

      _vector_add_assign(rs.s1, phi[ix]->s1);
      _vector_i_sub_assign(rs.s2, phi[ix]->s1);
    
      ix++;
      /*********************** direction -1 ************************/

      _su3_inverse_multiply(chi,(*U), phi[ix]->s0);
      _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);
      _complexcjg_times_vector(psi,ka1,chi);
      _complexcjg_times_vector(psi2,ka1,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_i_add_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_i_add_assign(rs.s2, psi2);

      U++;
      ix++;

      /*********************** direction +2 ************************/

      _vector_add_assign(rs.s0, phi[ix]->s0);
      _vector_add_assign(rs.s3, phi[ix]->s0);

      _vector_add_assign(rs.s1, phi[ix]->s1);
      _vector_sub_assign(rs.s2, phi[ix]->s1);
    
      ix++;
      /*********************** direction -2 ************************/

      _su3_inverse_multiply(chi,(*U), phi[ix]->s0);
      _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);
      _complexcjg_times_vector(psi,ka2,chi);
      _complexcjg_times_vector(psi2,ka2,chi2);
      _vector_add_assign(rs.s0, psi);
      _vector_sub_assign(rs.s3, psi);
      _vector_add_assign(rs.s1, psi2);
      _vector_add_assign(rs.s2, psi2);

      U++;
      ix++;
      /*********************** direction +3 ************************/

      _vector_add_assign(rs.s0, phi[ix]->s0);
      _vector_i_sub_assign(rs.s2, phi[ix]->s0);

      _vector_add_assign(rs.s1, phi[ix]->s1);
      _vector_i_add_assign(rs.s3, phi[ix]->s1);

      ix++;

      /*********************** direction -3 ************************/

      _su3_inverse_multiply(chi,(*U), phi[ix]->s0);
      _su3_inverse_multiply(chi2, (*U), phi[ix]->s1);
      _complexcjg_times_vector(psi,ka3,chi);
      _complexcjg_times_vector(psi2,ka3,chi2);      
      _vector_add(s->s0, rs.s0, psi);
      _vector_i_add(s->s2, rs.s2, psi);
      _vector_add(s->s1, rs.s1, psi2);
      _vector_i_sub(s->s3, rs.s3, psi2);

#ifndef OMP
      U++;
      ix++;
      s++;
#endif
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif

#ifdef OMP
  } /* OpenMP closing brace */
#endif

}
