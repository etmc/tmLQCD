/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "init_gauge_tmp.h"

su3 * gauge_tmp_ = NULL;
su3 ** gauge_tmp = NULL;

int init_gauge_tmp(const int V) {
  int i=0;

  gauge_tmp = calloc(V, sizeof(su3*));
  if(errno == ENOMEM) {
    return(1);
  }
  gauge_tmp_ = calloc(4*V+1, sizeof(su3));
  if(errno == ENOMEM) {
    return(1);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  gauge_tmp[0] = (su3*)(((unsigned long int)(gauge_tmp_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  gauge_tmp[0] = gauge_tmp_;
#endif
  for(i = 1; i < V; i++){
    gauge_tmp[i] = gauge_tmp[i-1]+4;
  }
  return(0);
}

void free_gauge_tmp() {
  free(gauge_tmp_);
  free(gauge_tmp);
}
