/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"

spinor * sp = NULL;

int init_spinor_field(const int V, const int nr) {
  int i = 0;

  sp = (spinor*)calloc(nr*V+1, sizeof(spinor));
  if(errno == ENOMEM) {
    return(1);
  }
  spinor_field = malloc(nr*sizeof(spinor*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  spinor_field[0] = (spinor*)(((unsigned long int)(sp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  spinor_field[0] = sp;
#endif
  
  for(i = 1; i < nr; i++){
    spinor_field[i] = spinor_field[i-1]+V;
  }

  return(0);
}

void free_spinor_field() {

  free(sp);
}
