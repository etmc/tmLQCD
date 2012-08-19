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
void Hopping_Matrix(int ieo, spinor * const l, spinor * const k){

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

#    if (defined MPI && !(defined _NO_COMM))
  xchange_field(k, ieo);
#    endif

#ifdef OMP
#pragma omp parallel
  {
#endif
  int ioff;
  int * hi;
  int ix,iy;
  int icy;
  su3 * restrict ALIGN up, * restrict ALIGN um;
  spinor * restrict rn, * restrict sp, * restrict sm;
  
#ifdef XLC
#pragma disjoint(temp, *sp, *sm, *rn, *up, *um)
#endif
  _declare_regs();

  if(ieo == 0){
    ioff = 0;
  } 
  else{
    ioff = (VOLUME+RAND)/2;
  } 
  /**************** loop over all lattice sites ****************/
#ifdef OMP
#pragma omp for
#endif
  for (int icx = ioff; icx < (VOLUME/2 + ioff); icx++){
    ix=g_eo2lexic[icx];

    rn=l+(icx-ioff);

    /*********************** direction +0 ************************/
    iy=g_iup[ix][0]; icy=g_lexic2eosub[iy];


    sp=k+icy;
#    if ((defined _GAUGE_COPY))
    up=&g_gauge_field_copy[icx][0];
#    else
    up=&g_gauge_field[ix][0];
#    endif
      
    _hop_t_p();

    /*********************** direction -0 ************************/

    iy=g_idn[ix][0]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    if ((defined _GAUGE_COPY))
    um = up+1;
#    else
    um=&g_gauge_field[iy][0];
#    endif

    _hop_t_m();

    /*********************** direction +1 ************************/

    iy=g_iup[ix][1]; icy=g_lexic2eosub[iy];

    sp=k+icy;

#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif

    _hop_x_p();      

    /*********************** direction -1 ************************/

    iy=g_idn[ix][1]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um=&g_gauge_field[iy][1];
#    else
    um=up+1;
#    endif

    _hop_x_m();

    /*********************** direction +2 ************************/

    iy=g_iup[ix][2]; icy=g_lexic2eosub[iy];

    sp=k+icy;
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif 

    _hop_y_p();

    /*********************** direction -2 ************************/

    iy=g_idn[ix][2]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][2];
#    else
    um = up +1;
#    endif

    _hop_y_m();

    /*********************** direction +3 ************************/

    iy=g_iup[ix][3]; icy=g_lexic2eosub[iy];

    sp=k+icy;
#    if ((defined _GAUGE_COPY))
    up=um+1;
#    else
    up+=1;
#    endif 

    _hop_z_p();

    /*********************** direction -3 ************************/

    iy=g_idn[ix][3]; icy=g_lexic2eosub[iy];

    sm=k+icy;
#    ifndef _GAUGE_COPY
    um = &g_gauge_field[iy][3];
#    else
    um = up+1;
#    endif

    _hop_z_m();

    _store_res();
  }
#ifdef OMP
  } /* OpenMP closing brace */
#endif

}
