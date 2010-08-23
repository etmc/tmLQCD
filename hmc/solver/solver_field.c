/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2009 Carsten Urbach
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
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<errno.h>
#include"global.h"
#include"su3.h"
#include"solver_field.h"

static int ini_sf = 0;
static int NR = 0;
static spinor * sf = NULL;

spinor ** solver_field = NULL;

int init_solver_field(const int V, const int nr){
  int i=0;

  if((ini_sf == 0) || (nr > NR)){
    ini_sf = 1;
    if(nr > NR){
      free(sf);
      free(solver_field);
      NR = nr;
    }
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
    if((void*)(sf = (spinor*)shmalloc((nr*V+1)*sizeof(spinor))) == NULL) {
      fprintf (stderr, "malloc errno in init_solver_field: %d\n",errno); 
      errno = 0;
      return(1);
    }
#else
    if((void*)(sf = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
      printf ("malloc errno in init_solver_field: %d\n",errno); 
      errno = 0;
      return(1);
    }
#endif
  if((void*)(solver_field = (spinor**)malloc(nr*sizeof(spinor*))) == NULL) {
    printf ("malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
    solver_field[0] = (spinor*)(((unsigned long int)(sf)+ALIGN_BASE)&~ALIGN_BASE);
#else
    solver_field[0] = sf;
#endif
    for(i = 1; i < nr; i++){
      solver_field[i] = solver_field[i-1]+V;
    }
  }
  return(0);
}

void finalize_solver(){
  if(ini_sf != 0){
    ini_sf = 0;
    free(sf);
    free(solver_field);
  }
  sf = NULL;
  solver_field = NULL;
}
