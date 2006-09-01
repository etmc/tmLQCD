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

spinor * sp = NULL;
spinor * sp_csg = NULL;

int init_spinor_field(const int V, const int nr) {
  int i = 0;

  sp = (spinor*)calloc(nr*V+1, sizeof(spinor));
  if(errno == ENOMEM) {
    return(1);
  }
  g_spinor_field = (spinor**)malloc(nr*sizeof(spinor*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_spinor_field[0] = (spinor*)(((unsigned long int)(sp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_spinor_field[0] = sp;
#endif
  
  for(i = 1; i < nr; i++){
    g_spinor_field[i] = g_spinor_field[i-1]+V;
  }

  return(0);
}

void free_spinor_field() {

  free(sp);
  free(sp_csg);
}


int init_csg_field(const int V, int * const nr) {
  int i = 0, j = 0;
  
  /* if all histories are zero, we do not need initialisation */
  if((nr[0] != 0) || (nr[2] != 0) || (nr[4] != 0) || (nr[6] != 0)) {
    sp_csg = (spinor*)calloc((nr[0]+nr[2]+nr[4]+nr[6])*V+1, sizeof(spinor));
    if(errno == ENOMEM) {
      return(1);
    }
    for(i = 0; i < 4; i++) {
      if(nr[2*i]!=0) {
	g_csg_field[i] = malloc(nr[2*i]*sizeof(spinor*));
	if(errno == ENOMEM) {
	  return(2);
	}
      }
      else g_csg_field[i] = NULL;
    }
#if ( defined SSE || defined SSE2 || defined SSE3)
    g_csg_field[0][0] = (spinor*)(((unsigned long int)(sp_csg)+ALIGN_BASE)&~ALIGN_BASE);
#else
    g_csg_field[0][0] = sp_csg;
#endif
    
    for(i = 1; i < nr[0]; i++){
      g_csg_field[0][i] = g_csg_field[0][i-1]+V;
    }
    for(j = 1; j < 4; j++) {
      if(nr[2*(j)]!=0) {
	g_csg_field[j][0] = g_csg_field[j-1][nr[2*(j-1)]-1]+V;
	for(i = 1; i < nr[2*j]; i++) {
	  g_csg_field[j][i] = g_csg_field[j][i-1]+V;
	}
      }
    }
    
    g_csg_index_array[0] = (int*) malloc((nr[0]+nr[2]+nr[4]+nr[6])*sizeof(int));
    for(i = 1; i < 4; i++) {
      g_csg_index_array[i] = g_csg_index_array[i-1]+nr[2*(i-1)];
    }
  }

  return(0);
}
