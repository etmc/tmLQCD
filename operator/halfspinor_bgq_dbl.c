/**********************************************************************
 *
 *
 * Copyright (C) 2012 Carsten Urbach
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
  int ix;
  su3 * restrict ALIGN U;
  spinor * restrict ALIGN s;
  halfspinor * restrict * phi ALIGN;
  halfspinor32 * restrict * phi32 ALIGN;
  /* We have 32 registers available */
  _declare_hregs();

#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif
#pragma disjoint(*s, *U)

#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

  __alignx(16, l);
  __alignx(16, k);
  if(g_sloppy_precision == 1 && g_sloppy_precision_flag == 1) {
    __alignx(16, HalfSpinor32);
    /* We will run through the source vector now */
    /* instead of the solution vector            */
    s = k;
    _prefetch_spinor(s); 

    /* s contains the source vector */

    if(ieo == 0) {
      U = g_gauge_field_copy[0][0];
    }
    else {
      U = g_gauge_field_copy[1][0];
    }
    phi32 = NBPointer32[ieo];

    _prefetch_su3(U);
    /**************** loop over all lattice sites ******************/
    ix=0;
    for(int i = 0; i < (VOLUME)/2; i++){
      /*********************** direction +0 ************************/
      _hop_t_p_pre32();
      s++; 
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _hop_t_m_pre32();
      ix++;

      /*********************** direction +1 ************************/
      _hop_x_p_pre32();
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _hop_x_m_pre32();
      ix++;

      /*********************** direction +2 ************************/
      _hop_y_p_pre32();

      ix++;
      U++;

      /*********************** direction -2 ************************/
      _hop_y_m_pre32();
      ix++;

      /*********************** direction +3 ************************/
      _hop_z_p_pre32();
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _hop_z_m_pre32();
      ix++;

      /************************ end of loop ************************/
    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield32(); 
#    endif
    s = l;
    phi32 = NBPointer32[2 + ieo];
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
    //_prefetch_halfspinor(phi32[0]);
    _prefetch_su3(U);
  
    /* Now we sum up and expand to a full spinor */
    ix = 0;
    /*   _prefetch_spinor_for_store(s); */
    for(int i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      /*     _prefetch_spinor_for_store(s); */
      //_prefetch_halfspinor(phi32[ix+1]);
      /*********************** direction +0 ************************/
      _hop_t_p_post32();
      ix++;
      /*********************** direction -0 ************************/
      _hop_t_m_post32();
      U++;
      ix++;
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
      U++;
      ix++;
      s++;
    }
  }
  else {
    __alignx(16, HalfSpinor);
    /* We will run through the source vector now */
    /* instead of the solution vector            */
    s = k;
    _prefetch_spinor(s); 

    /* s contains the source vector */

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
    for(int i = 0; i < (VOLUME)/2; i++){
      /*********************** direction +0 ************************/
      _hop_t_p_pre();
      s++; 
      U++;
      ix++;

      /*********************** direction -0 ************************/
      _hop_t_m_pre();
      ix++;

      /*********************** direction +1 ************************/
      _hop_x_p_pre();
      ix++;
      U++;

      /*********************** direction -1 ************************/
      _hop_x_m_pre();
      ix++;


      /*********************** direction +2 ************************/
      _hop_y_p_pre();
      ix++;
      U++;

      /*********************** direction -2 ************************/
      _hop_y_m_pre();
      ix++;

      /*********************** direction +3 ************************/
      _hop_z_p_pre();
      ix++;
      U++;

      /*********************** direction -3 ************************/
      _hop_z_m_pre();
      ix++;

      /************************ end of loop ************************/

    }

#    if (defined MPI && !defined _NO_COMM)
    xchange_halffield(); 
#    endif
    s = l;
    phi = NBPointer[2 + ieo];
    //_prefetch_halfspinor(phi[0]);
    if(ieo == 0) {
      U = g_gauge_field_copy[1][0];
    }
    else {
      U = g_gauge_field_copy[0][0];
    }
    _prefetch_su3(U);
  
    /* Now we sum up and expand to a full spinor */
    ix = 0;
    /*   _prefetch_spinor_for_store(s); */
    for(int i = 0; i < (VOLUME)/2; i++){
      /* This causes a lot of trouble, do we understand this? */
      /*     _prefetch_spinor_for_store(s); */
      //_prefetch_halfspinor(phi[ix+1]);
      /*********************** direction +0 ************************/
      _hop_t_p_post();
      ix++;
      /*********************** direction -0 ************************/
      _hop_t_m_post();
      U++;
      ix++;
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
      U++;
      ix++;
      s++;
    }
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif
}
