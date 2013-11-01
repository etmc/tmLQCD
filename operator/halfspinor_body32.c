/**********************************************************************
 * single precision version Copyright (C) 2013 Florian Burger
 * based on halfspinor_body.c by Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator 
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002 
 * this is a new version based on the aforementioned implementations
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


int ix;
su3_32 * restrict U ALIGN32;
spinor32 * restrict s ALIGN32;
halfspinor32 * restrict * phi2 ALIGN32;
_declare_hregs();

#ifdef XLC
# pragma disjoint(*l, *k)
# pragma disjoint(*k, *U)
# pragma disjoint(*l, *U)
# pragma disjoint(*U, *s)
# pragma disjoint(*k, *s)
# pragma disjoint(*l, *s)
__alignx(16, l);
__alignx(16, k);
__alignx(16, U);
__alignx(16, s);
#endif 

//convert kappas to float locally
_Complex float ALIGN32 ka0_32 = (_Complex float) ka0;
_Complex float ALIGN32 ka1_32 = (_Complex float) ka1;
_Complex float ALIGN32 ka2_32 = (_Complex float) ka2;
_Complex float ALIGN32 ka3_32 = (_Complex float) ka3;

#ifndef OMP  
s = k;
_prefetch_spinor_32(s);
 U = g_gauge_field_copy_32[ieo][0];
_prefetch_su3_32(U);
#else
  u0 = g_gauge_field_copy_32[ieo][0];
#endif

  phi2 = NBPointer32[ieo];
  
#ifdef OMP
#pragma omp for
#else
  ix=0;
#endif
  for(unsigned int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    s=k+i;
    _prefetch_spinor_32(s);    
    ix=i*8;
    U=u0+i*4;    
    _prefetch_su3_32(U);    
#endif
    _hop_t_p_pre32();
    U++;
    ix++;
    
    _hop_t_m_pre32();
    ix++;
    
    _hop_x_p_pre32();
    U++;
    ix++;
    
    _hop_x_m_pre32();
    ix++;
    
    _hop_y_p_pre32();
    U++;
    ix++;
    
    _hop_y_m_pre32();
    ix++;
    
    _hop_z_p_pre32();
    U++;
    ix++;
    
    _hop_z_m_pre32();
    
#ifndef OMP
    s++;
    ix++;
#endif
  }

#if (defined MPI && !defined _NO_COMM) 
#  ifdef OMP
#  pragma omp single nowait
   {
#  endif
 
#  ifdef SPI

     // Initialize the barrier, resetting the hardware.
     int rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route*/  );
     if(rc) {
       printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
       exit(__LINE__);
     }
     // reset the recv counter 
     recvCounter = totalMessageSize/2;
     global_barrier(); // make sure everybody is set recv counter

     //#pragma omp for nowait
     for (unsigned int j = 0; j < spi_num_dirs; j++) {
       descCount[ j ] =
	 msg_InjFifoInject ( injFifoHandle,
			     j,
			     &SPIDescriptors32[j]);
     }

#  else // SPI
     xchange_halffield32(); 
#    ifdef _THREAD_OVERLAP     
        wait_halffield32();     
#    endif

#  endif  // SPI
    
#  ifdef OMP
   }
#  endif

#endif // MPI
 
#ifndef OMP
  s = l;
  if(ieo == 0) {
    U = g_gauge_field_copy_32[1][0];
  }
  else {
    U = g_gauge_field_copy_32[0][0];
  }
 _prefetch_su3_32(U);
#else
  if(ieo == 0) {
    u0 = g_gauge_field_copy_32[1][0];
  }
  else {
    u0 = g_gauge_field_copy_32[0][0];
  }
#endif
  
  phi2 = NBPointer32[2 + ieo];
  
#ifdef OMP
#ifdef _THREAD_OVERLAP  
#pragma omp for schedule(guided,32)
#else
#pragma omp for
#endif
#else
  ix = 0;
#endif
  for(unsigned int i = 0; i < bodyV; i++){
#ifdef OMP
    ix=i*8;
    U=u0+i*4;  
    _prefetch_su3_32(U);    
    s=l+i;
    _prefetch_spinor_32(s);
#endif
#ifdef _TM_SUB_HOP
     pn=p+i;
#endif
    _prefetch_halfspinor(phi2[ix+1]);     
    _hop_t_p_post32();
    ix++;
    
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_t_m_post32();
    ix++;
    U++;
    
    _prefetch_halfspinor(phi2[ix+1]);
    _hop_x_p_post32();
    ix++;
    
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_x_m_post32();
    U++;
    ix++;
    
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_y_p_post32();
    ix++;
    
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_y_m_post32();
    U++;
    ix++;
   
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_z_p_post32();
    ix++;
   
    _prefetch_halfspinor(phi2[ix+1]);    
    _hop_z_m_post32();
    
#ifdef _MUL_G5_CMPLX
    _hop_mul_g5_cmplx_and_store32(s);
#elif defined _TM_SUB_HOP
     _g5_cmplx_sub_hop_and_g5store32(s);
#else
    _hop_store_post32(s);
#endif
    
#ifndef OMP
    U++;
    ix++;
    s++;
#endif
  }
 
#ifndef _THREAD_OVERLAP
#ifdef OMP
#pragma omp single
   {
#endif
     
#    if (defined MPI && !defined _NO_COMM)
#      ifdef SPI
     // wait for receive completion
     while ( recvCounter > 0 );
     _bgq_msync();
#      else
     wait_halffield32();
#      endif
#    endif
     
#ifdef OMP
   }
#endif
#endif /* _THREAD_OVERLAP */


#ifdef OMP
#pragma omp for
#endif
   /* #pragma ivdep */
   for(unsigned int i = bodyV; i < VOLUME/2; i++){
#ifdef OMP
     ix=i*8;
     U=u0+i*4;
     _prefetch_su3_32(U);
     s=l+i;
     _prefetch_spinor_32(s);
#endif
#ifdef _TM_SUB_HOP
     pn=p+i;
#endif
     _prefetch_halfspinor(phi2[ix+1]);    
     _hop_t_p_post32();
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);     
     _hop_t_m_post32();
     ix++;
     U++;
     _prefetch_halfspinor(phi2[ix+1]);
     _hop_x_p_post32();
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);          
     _hop_x_m_post32();
     U++;
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);     
     _hop_y_p_post32();
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);     
     _hop_y_m_post32();
     U++;
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);     
     _hop_z_p_post32();
     ix++;
     _prefetch_halfspinor(phi2[ix+1]);     
     _hop_z_m_post32();

#ifdef _MUL_G5_CMPLX
     _hop_mul_g5_cmplx_and_store32(s);
#elif defined _TM_SUB_HOP
     _g5_cmplx_sub_hop_and_g5store32(s);
#else
     _hop_store_post32(s);
#endif
     
#ifndef OMP
     U++;
     ix++;
     s++;
#endif
   }
 

