/**********************************************************************
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



//  if(k == l){
//    printf("Error in H_psi (simple.c):\n");
//    printf("Arguments k and l must be different\n");
//    printf("Program aborted\n");
//    exit(1);
//  }


int ix;
su3 * restrict U ALIGN;
spinor * restrict s ALIGN;
halfspinor * restrict * phi ALIGN;
halfspinor32 * restrict * phi32 ALIGN;
_declare_hregs();

#ifdef XLC
# pragma disjoint(*l, *k)
# pragma disjoint(*k, *U)
# pragma disjoint(*l, *U)
# pragma disjoint(*U, *s)
# pragma disjoint(*k, *s)
# pragma disjoint(*l, *s)
__alignx(32, l);
__alignx(32, k);
__alignx(32, U);
__alignx(32, s);
#endif


#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif

#ifndef OMP  
s = k;
_prefetch_spinor(s);
if(ieo == 0) {
  U = g_gauge_field_copy[0][0];
 }
 else {
   U = g_gauge_field_copy[1][0];
 }
_prefetch_su3(U);
#else
if(ieo == 0) {
  u0 = g_gauge_field_copy[0][0];
 }
 else {
   u0 = g_gauge_field_copy[1][0];
 }
#endif
#if (defined SSE2 || defined SSE3)
g_sloppy_precision = 0;
#endif
if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
  phi32 = NBPointer32[ieo];
  
  /**************** loop over all lattice sites ****************/
#ifdef OMP
#pragma omp for
#else
  ix=0;
#endif
  for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    U=u0+i*4;
    s=k+i;
    ix=i*8;
#endif
    /*********************** direction +0 ************************/
    _hop_t_p_pre32();
    U++;
    ix++;
    
    /*********************** direction -0 ************************/
    _hop_t_m_pre32();
    ix++;
    
    /*********************** direction +1 ************************/
    _hop_x_p_pre32();
    U++;
    ix++;
    
    /*********************** direction -1 ************************/
    _hop_x_m_pre32();
    ix++;
    
    /*********************** direction +2 ************************/
    _hop_y_p_pre32();
    U++;
    ix++;
    
    /*********************** direction -2 ************************/
    _hop_y_m_pre32();
    ix++;
    
    /*********************** direction +3 ************************/
    _hop_z_p_pre32();
    U++;
    ix++;
    
    /*********************** direction -3 ************************/
    _hop_z_m_pre32();
    
#ifndef OMP
    s++;
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
    u0 = g_gauge_field_copy[1][0];
  }
  else {
    u0 = g_gauge_field_copy[0][0];
  }
#endif
  
  phi32 = NBPointer32[2 + ieo];
  
#ifdef OMP
#pragma omp for
#else
  ix = 0;
#endif
  for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    ix=i*8;
    s=l+i;
    U=u0+i*4;
#endif
    /*********************** direction +0 ************************/
    _hop_t_p_post32();
    ix++;
    
    /*********************** direction -0 ************************/
    _hop_t_m_post32();
    ix++;
    U++;
    
    /*********************** direction +1 ************************/
    _hop_x_p_post32();
    ix++;
    
    /*********************** direction -1 ************************/
    _hop_x_m_post32();
    U++;
    ix++;
    
    /*********************** direction +2 ************************/
    _hop_y_p_post32();
    ix++;
    
    /*********************** direction -2 ************************/
    _hop_y_m_post32();
    U++;
    ix++;
    
    /*********************** direction +3 ************************/
    _hop_z_p_post32();
    ix++;
    
    /*********************** direction -3 ************************/
    _hop_z_m_post32();
    
#ifdef _MUL_G5_CMPLX
    _hop_mul_g5_cmplx_and_store(s);
#else
    _hop_store_post(s);
#endif
    
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
   for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
     s=k+i;
     _prefetch_spinor(s);
     ix=i*8;
     U=u0+i*4;
     _prefetch_su3(U);
#endif
     
     /*********************** direction +0 ************************/
     _hop_t_p_pre();
     U++;
     ix++;
     
     /*********************** direction -0 ************************/
     _hop_t_m_pre();
     ix++;
     
     /*********************** direction +1 ************************/
     _hop_x_p_pre();
     U++;
     ix++;
     
     /*********************** direction -1 ************************/
     _hop_x_m_pre();
     ix++;
     
     /*********************** direction +2 ************************/
     _hop_y_p_pre();
     U++;
     ix++;
     
     /*********************** direction -2 ************************/
     _hop_y_m_pre();
     ix++;
     
     /*********************** direction +3 ************************/
     _hop_z_p_pre();
     U++;
     ix++;
     
     /*********************** direction -3 ************************/
     _hop_z_m_pre();
     
#ifndef OMP
     s++;            
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
   _prefetch_su3(U);
#else
   if(ieo == 0) {
     u0 = g_gauge_field_copy[1][0];
   }
   else {
     u0 = g_gauge_field_copy[0][0];
   }
#endif
   
   phi = NBPointer[2 + ieo];
   
#ifdef OMP
#pragma omp for
#else
   ix = 0;
#endif
   /* #pragma ivdep */
   for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
     ix=i*8;
     U=u0+i*4;
     _prefetch_su3(U);
     s=l+i;
     _prefetch_spinor(s);
#endif
     
     /*********************** direction +0 ************************/
     _hop_t_p_post();
     ix++;
     
     /*********************** direction -0 ************************/
     _hop_t_m_post();
     ix++;
     U++;
     
     /*********************** direction +1 ************************/
     _hop_x_p_post();
     ix++;
     
     /*********************** direction -1 ************************/
     _hop_x_m_post();
     U++;
     ix++;
     
     /*********************** direction +2 ************************/
     _hop_y_p_post();
     ix++;
     
     /*********************** direction -2 ************************/
     _hop_y_m_post();
     U++;
     ix++;
     
     /*********************** direction +3 ************************/
     _hop_z_p_post();
     ix++;
     
     /*********************** direction -3 ************************/
     _hop_z_m_post();
     
#ifdef _MUL_G5_CMPLX
     _hop_mul_g5_cmplx_and_store(s);
#else
     _hop_store_post(s);
#endif
     
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



