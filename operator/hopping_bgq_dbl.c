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

void Hopping_Matrix(const int ieo, spinor * const l, spinor * const k) {
#ifdef _GAUGE_COPY
  if(g_update_gauge_copy) {
    update_backward_gauge(g_gauge_field);
  }
#endif

#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

#ifdef OMP
#pragma omp parallel
  {
#endif
  int icx,icy,icz,ioff,ioff2;
  int ix,iy,iz;
  su3 * restrict ALIGN up;
  su3 * restrict ALIGN um;
  spinor * restrict ALIGN sp;
  spinor * restrict ALIGN sm;
  spinor * restrict ALIGN rn;
  
#pragma disjoint(*sp, *sm, *rn, *up, *um, *l, *k)
  _declare_regs();

  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  }

#ifndef OMP
  ix=g_eo2lexic[ioff];
  iy=g_iup[ix][0]; 
  icy=g_lexic2eosub[iy];

  sp=k+icy;

#    if ((defined _GAUGE_COPY))
  up=&g_gauge_field_copy[ioff][0];
#    else
  up=&g_gauge_field[ix][0];
#    endif
#endif

  /**************** loop over all lattice sites ******************/
#ifdef OMP
#pragma omp for
#endif
  for(icx = ioff; icx < (VOLUME/2+ioff); icx++){

#ifdef OMP
    ix=g_eo2lexic[icx];
    iy=g_iup[ix][0]; 
    icy=g_lexic2eosub[iy];

    sp=k+icy;
#endif
#    if ((defined _GAUGE_COPY))
    p=&g_gauge_field_copy[icx][0];
#    else
    p=&g_gauge_field[ix][0];
#    endif
#endif

    rn=l+(icx-ioff);
#ifndef OMP
    ix=g_eo2lexic[icx];
#endif
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

    _store_res();
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif
  return;
}
