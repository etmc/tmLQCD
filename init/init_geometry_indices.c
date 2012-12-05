/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "init_geometry_indices.h"

int *iup = NULL, *idn = NULL, *ipt = NULL, **ipt_ = NULL, ***ipt__ = NULL;

int init_geometry_indices(const int V) {
  int i = 0;

  g_idn= (int**)calloc(V, sizeof(int*));
  if((void*)g_idn == NULL) return(1);
  g_iup = (int**)calloc(V, sizeof(int*));
  if((void*)g_iup == NULL) return(2);

  idn = (int*)calloc(4*V, sizeof(int));
  if((void*)idn == NULL ) return(6);
  iup = (int*)calloc(4*V, sizeof(int));
  if((void*)iup == NULL) return(7);

  g_ipt = (int****)calloc(T+4,sizeof(int*));
  if((void*)g_ipt == NULL) return(5);
  ipt__ = (int***)calloc ((T+4)*(LX+4), sizeof(int*));
  if((void*)ipt__ == NULL) return(4);
  ipt_ = (int**)calloc((T+4)*(LX+4)*(LY+4), sizeof(int*));
  if((void*)ipt_ == NULL) return(3);
  ipt = (int*)calloc((T+4)*(LX+4)*(LY+4)*(LZ+4), sizeof(int));
  if((void*)ipt == NULL) return(8);

  g_lexic2eo = (int*)calloc(V, sizeof(int));
  if((void*)g_lexic2eo == NULL) return(9);
  /* this +2 is for sanity reasons */
  g_lexic2eosub = (int*)calloc(V+2, sizeof(int));
  if((void*)g_lexic2eosub == NULL) return(10);
  g_eo2lexic = (int*)calloc(V, sizeof(int));
  if((void*)g_eo2lexic == NULL) return(11);

#if ( defined PARALLELXYZT || defined PARALLELXYZ )
  g_field_z_ipt_even = (int*)calloc(T*LX*LY, sizeof(int));
  if((void*)g_field_z_ipt_even == NULL) return(12);
  g_field_z_ipt_odd  = (int*)calloc(T*LX*LY, sizeof(int));
  if((void*)g_field_z_ipt_odd == NULL) return(13);

  g_field_z_disp_even_dn = (int*)calloc(T*LX*LY/2, sizeof(int));
  if((void*)g_field_z_disp_even_dn == NULL) return(14);
  g_field_z_disp_even_up = (int*)calloc(T*LX*LY/2, sizeof(int));
  if((void*)g_field_z_disp_even_up == NULL) return(15);
  g_field_z_disp_odd_dn = (int*)calloc(T*LX*LY/2, sizeof(int));
  if((void*)g_field_z_disp_odd_dn == NULL) return(16);
  g_field_z_disp_odd_up = (int*)calloc(T*LX*LY/2, sizeof(int));
  if((void*)g_field_z_disp_odd_up == NULL) return(17);
#endif

#ifdef _USE_TSPLITPAR
  g_1st_eot= (int**)calloc(T, sizeof(int*));
  if((void*)g_1st_eot == NULL) return(18);
  for(i=0;i<T;i++){
    g_1st_eot[i]= (int *)calloc(2, sizeof(int));
  }
  g_1st_xt_int_dn = (int*)calloc(T, sizeof(int));
  g_1st_xt_int_up = (int*)calloc(T, sizeof(int));
  g_1st_xt_ext_dn = (int*)calloc(T, sizeof(int));
  g_1st_xt_ext_up = (int*)calloc(T, sizeof(int));
  g_1st_yt_int_dn = (int*)calloc(T, sizeof(int));
  g_1st_yt_int_up = (int*)calloc(T, sizeof(int));
  g_1st_yt_ext_dn = (int*)calloc(T, sizeof(int));
  g_1st_yt_ext_up = (int*)calloc(T, sizeof(int));
  g_1st_zt_int_dn = (int*)calloc(T, sizeof(int));
  g_1st_zt_int_up = (int*)calloc(T, sizeof(int));
  g_1st_zt_ext_dn = (int*)calloc(T, sizeof(int));
  g_1st_zt_ext_up = (int*)calloc(T, sizeof(int));

  g_field_zt_disp_even_dn = (int **)calloc(T, sizeof(int*));
  g_field_zt_disp_even_up = (int **)calloc(T, sizeof(int*));
  g_field_zt_disp_odd_dn = (int **)calloc(T, sizeof(int*));
  g_field_zt_disp_odd_up = (int **)calloc(T, sizeof(int*));
  for(i=0;i<T;i++){
    g_field_zt_disp_even_dn[i] = (int *)calloc((LX*LY+1)/2, sizeof(int));
    g_field_zt_disp_even_up[i] = (int *)calloc((LX*LY+1)/2, sizeof(int));
    g_field_zt_disp_odd_dn[i] = (int *)calloc((LX*LY+1)/2, sizeof(int));
    g_field_zt_disp_odd_up[i] = (int *)calloc((LX*LY+1)/2, sizeof(int));
  }

#endif

  g_coord= (int**)calloc(VOLUME, sizeof(int*));
  if((void*)g_coord == NULL) return(19);
  for(i=0;i<VOLUME;i++){
    g_coord[i]= (int*)calloc(4, sizeof(int));
  }

  g_iup_eo= (int**)calloc(VOLUME+RAND, sizeof(int*));  // NEW GIUPDNEO
  if((void*)g_iup_eo == NULL) return(21);
  for(i=0;i<VOLUME+RAND;i++){
    g_iup_eo[i]= (int*)calloc(4, sizeof(int));
  }

  g_idn_eo= (int**)calloc(VOLUME+RAND, sizeof(int*));
  if((void*)g_idn_eo == NULL) return(22);
  for(i=0;i<VOLUME+RAND;i++){
    g_idn_eo[i]= (int*)calloc(4, sizeof(int));
  }


  /* This should only be used for the SFBC. */
  /* This should not be used for anything other than the SFBC */
  /* because it might eventually vanish. */

  g_idn[0] = idn;
  g_iup[0] = iup;

  ipt_[0] = ipt;
  ipt__[0] = ipt_;
  g_ipt[0] = ipt__;
  for(i = 1; i < V; i++){
    g_idn[i] = g_idn[i-1]+4;
    g_iup[i] = g_iup[i-1]+4;
  }
  for(i = 1; i < (T+4)*(LX+4)*(LY+4); i++){
    ipt_[i] = ipt_[i-1]+(LZ+4);
  }
  for(i = 1; i < (T+4)*(LX+4); i++){
    ipt__[i] = ipt__[i-1]+(LY+4);
  }
  for(i = 1; i < (T+4); i++){
    g_ipt[i] = g_ipt[i-1]+(LX+4);
  }
  g_hi = (int*)calloc(16*(VOLUME+RAND)+2,sizeof(int));
  if((void*) g_hi == NULL) return(40);

#ifdef WITHLAPH
  g_idn3d = (int**)calloc(SPACEVOLUME, sizeof(int*));
  if((void*)g_idn == NULL) return(31);
  g_iup3d = (int**)calloc(SPACEVOLUME, sizeof(int*));
  if((void*)g_iup == NULL) return(32);
  for (i=0;i<SPACEVOLUME;i++){
    g_idn3d[i] = (int*)calloc(4, sizeof(int));
    g_iup3d[i] = (int*)calloc(4, sizeof(int));
  }
#endif

  return(0);
}

void free_geometry_indices() {
  free(idn); 
  free(iup); 
  free(ipt);
  free(ipt_);
  free(ipt__);
  free(g_ipt);
  free(g_hi);
  free(g_idn);
  free(g_iup);
  free(g_eo2lexic);
  free(g_lexic2eosub);
  free(g_lexic2eo);
#if ( defined PARALLELXYZT || defined PARALLELXYZ )
  free(g_field_z_ipt_odd);
  free(g_field_z_ipt_even);
#endif
}
