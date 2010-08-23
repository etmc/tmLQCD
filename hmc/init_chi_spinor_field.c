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
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "init_chi_spinor_field.h"

spinor * sp_up = NULL;
spinor * sp_dn = NULL;

static int chi_initialised = 0;

int init_chi_spinor_field(const int V, const int nr) {
  int i = 0;
  static int _nr = 0;

  if(!chi_initialised || nr > _nr) {
    free_chi_spinor_field();
    _nr = nr;
    chi_initialised = 1;
    if((void*)(sp_up = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(g_chi_up_spinor_field = malloc(nr*sizeof(spinor*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(2);
    }
    if((void*)(sp_dn = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(g_chi_dn_spinor_field = malloc(nr*sizeof(spinor*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(2);
    }
#if ( defined SSE || defined SSE2 || defined SSE3)
    g_chi_up_spinor_field[0] = (spinor*)(((unsigned long int)(sp_up)+ALIGN_BASE)&~ALIGN_BASE);
    g_chi_dn_spinor_field[0] = (spinor*)(((unsigned long int)(sp_dn)+ALIGN_BASE)&~ALIGN_BASE);
#else
    g_chi_up_spinor_field[0] = sp_up;
    g_chi_dn_spinor_field[0] = sp_dn;
#endif
    
    for(i = 1; i < nr; i++){
      g_chi_up_spinor_field[i] = g_chi_up_spinor_field[i-1]+V;
      g_chi_dn_spinor_field[i] = g_chi_dn_spinor_field[i-1]+V;
    }
  }
  return(0);
}

void free_chi_spinor_field() {
  if(chi_initialised) {
    free(sp_up);
    free(sp_dn);
    free(g_chi_dn_spinor_field);
    free(g_chi_up_spinor_field);
  }
  return;
}


