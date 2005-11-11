/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"


bispinor * bisp = NULL;

int init_bispinor_field(const int V, const int nr) {
  int i = 0;

  bisp = (bispinor*)calloc(nr*V+1, sizeof(bispinor));
  if(errno == ENOMEM) {
    return(1);
  }
  g_bispinor_field = malloc(nr*sizeof(bispinor*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_bispinor_field[0] = (bispinor*)(((unsigned long int)(bisp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_bispinor_field[0] = bisp;
#endif
  
  for(i = 1; i < nr; i++){
    g_bispinor_field[i] = g_bispinor_field[i-1]+V;
  }

  return(0);
}

void free_bispinor_field() {

  free(bisp);
  /*  free(sp_csg); */
}

