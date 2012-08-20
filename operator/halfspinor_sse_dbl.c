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

/* input on k; output on l */
void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k){
#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

#ifdef OMP
#pragma omp parallel
{
  su3 * restrict U0 ALIGN;
#endif

  int ix;
  su3 * restrict U ALIGN;
  spinor * restrict s ALIGN;
  halfspinor ** phi ALIGN;
  _declare_hregs();

#ifdef _KOJAK_INST
#pragma pomp inst begin(hoppingmatrix)
#endif

#ifndef OMP
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
  _prefetch_su3(U);
#else
  if(ieo == 0) {
    U0 = g_gauge_field_copy[0][0];
  }
  else {
    U0 = g_gauge_field_copy[1][0];
  }
#endif
  phi = NBPointer[ieo];

  /**************** loop over all lattice sites ******************/
#ifdef OMP
#pragma omp for
#else
  ix = 0;
#endif
  for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    s = k+i;
    _prefetch_spinor(s);
    U = U0+i*4;
    _prefetch_su3(U);
    ix = i*8;
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
#ifndef OMP
    ix++;
    s++;
#endif
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
    U0 = g_gauge_field_copy[1][0];
  }
  else {
    U0 = g_gauge_field_copy[0][0];
  }
#endif
  phi = NBPointer[2 + ieo];

  
  /* Now we sum up and expand to a full spinor */
#ifdef OMP
#pragma omp for
#else
  ix = 0;
#endif
  for(int i = 0; i < (VOLUME)/2; i++){
#ifdef OMP
    U = U0 + i*4;
    _prefetch_su3(U);
    ix = i*8;
    s = l + i;
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
    ix++;
    U++;

    /*********************** direction +2 ************************/
    _hop_y_p_post();
    ix++;

    /*********************** direction -2 ************************/
    _hop_y_m_post();
    ix++;
    U++;
    /*********************** direction +3 ************************/
    _hop_z_p_post();

    ix++;
    /*********************** direction -3 ************************/
    _hop_z_m_post();
#ifndef OMP
    ix++;
    U++;
    s++;
#endif
  }
#ifdef _KOJAK_INST
#pragma pomp inst end(hoppingmatrix)
#endif

#ifdef OMP
  } /* omp parallel closing bracket */
#endif
}


