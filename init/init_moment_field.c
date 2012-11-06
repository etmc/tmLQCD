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
#include "su3.h"
#include "su3adj.h"
#include "sse.h"

#include <buffers/adjoint.h>


adjoint_field_t df, du, mo;

int init_moment_field(const int V, const int VR) {
  int i = 0;

  df = get_adjoint_field();
  du = get_adjoint_field();
  mo = get_adjoint_field();
  
  if((void*)(moment = (su3adj**)calloc(V,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
  
  moment[0] = (su3adj*)mo[0];
  for(i = 1; i < V; i++)
    moment[i] = moment[i-1]+4;
 
  if((void*)(df0 = (su3adj**)calloc(VR,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(4);
  }
  
  df0[0] = (su3adj*)df[0];
  for(i = 1; i < VR; i++)
    df0[i] = df0[i-1]+4; 
  
  if((void*)(ddummy = (su3adj**)calloc(VR,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(6);
  }

  ddummy[0] = (su3adj*)du[0];
  for(i = 1; i < VR; i++){
    ddummy[i] = ddummy[i-1]+4;
  }

  return(0);
}

void free_moment_field()
{
  return_adjoint_field(&df);
  return_adjoint_field(&du);
  return_adjoint_field(&mo);
  
  free(df0);
  free(ddummy);
  free(moment);
}
