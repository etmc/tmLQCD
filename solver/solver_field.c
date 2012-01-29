/***********************************************************************
 *
 * Copyright (C) 2009,2011 Carsten Urbach
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

int init_solver_field(spinor *** const solver_field, const int V, const int nr) {
  int i=0;

  /* allocate nr+1 to save the linear field in solver_field[nr] */
  if((void*)((*solver_field) = (spinor**)malloc((nr+1)*sizeof(spinor*))) == NULL) {
    printf ("malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(2);
  }
  
  /* allocate the full chunk of memory to solver_field[nr] */
#if (defined _USE_SHMEM && !(defined _USE_HALFSPINOR))
  if((void*)((*solver_field)[nr] = (spinor*)shmalloc((nr*V+1)*sizeof(spinor))) == NULL) {
    fprintf (stderr, "malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(1);
  }
#else
  if((void*)((*solver_field)[nr] = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(1);
  }
#endif

  /* now cut in pieces and distribute to solver_field[0]-solver_field[nr-1] */
#if ( defined SSE || defined SSE2 || defined SSE3)
  (*solver_field)[0] = (spinor*)(((unsigned long int)((*solver_field)[nr])+ALIGN_BASE)&~ALIGN_BASE);
#else
  (*solver_field)[0] = (*solver_field)[nr];
#endif
  for(i = 1; i < nr; i++){
    (*solver_field)[i] = (*solver_field)[i-1]+V;
  }
  return(0);
}

void finalize_solver(spinor ** solver_field, const int nr){
  free(solver_field[nr]);
  free(solver_field);
  solver_field = NULL;
}


int init_bisolver_field(bispinor *** const solver_field, const int V, const int nr) {
  int i=0;

  /* allocate nr+1 to save the linear field in solver_field[nr] */
  if((void*)((*solver_field) = (bispinor**)malloc((nr+1)*sizeof(bispinor*))) == NULL) {
    printf ("malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(2);
  }
  
  /* allocate the full chunk of memory to solver_field[nr] */
  if((void*)((*solver_field)[nr] = (bispinor*)calloc(nr*V+1, sizeof(bispinor))) == NULL) {
    printf ("malloc errno in init_solver_field: %d\n",errno); 
    errno = 0;
    return(1);
  }

  /* now cut in pieces and distribute to solver_field[0]-solver_field[nr-1] */
#if ( defined SSE || defined SSE2 || defined SSE3)
  (*solver_field)[0] = (bispinor*)(((unsigned long int)((*solver_field)[nr])+ALIGN_BASE)&~ALIGN_BASE);
#else
  (*solver_field)[0] = (*solver_field)[nr];
#endif
  for(i = 1; i < nr; i++){
    (*solver_field)[i] = (*solver_field)[i-1]+V;
  }
  return(0);
}

void finalize_bisolver(bispinor ** solver_field, const int nr) {
  free(solver_field[nr]);
  free(solver_field);
  solver_field = NULL;
}
