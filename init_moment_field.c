/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "sse.h"

su3adj * mo=NULL, *df=NULL, *du=NULL;

int init_moment_field(const int V, const int VR) {
  int i = 0;

  mo = (su3adj*)calloc(4*V+1, sizeof(su3adj));
  if(errno == ENOMEM) {
    return(1);
  }
  moment = malloc(V*sizeof(su3adj*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  moment[0] = (su3adj*)(((unsigned long int)(mo)+ALIGN_BASE)&~ALIGN_BASE);
#else
  moment[0] = mo;
#endif
  
  for(i = 1; i < V; i++){
    moment[i] = moment[i-1]+4;
  }

  df = (su3adj*)calloc(4*VR+1, sizeof(su3adj));
  if(errno == ENOMEM) {
    return(1);
  }
  df0 = malloc(VR*sizeof(su3adj*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  df0[0] = (su3adj*)(((unsigned long int)(df)+ALIGN_BASE)&~ALIGN_BASE);
#else
  df0[0] = df;
#endif
  
  for(i = 1; i < VR; i++){
    df0[i] = df0[i-1]+4;
  }

  du = (su3adj*)calloc(4*VR+1, sizeof(su3adj));
  if(errno == ENOMEM) {
    return(1);
  }
  ddummy = malloc(VR*sizeof(su3adj*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  ddummy[0] = (su3adj*)(((unsigned long int)(du)+ALIGN_BASE)&~ALIGN_BASE);
#else
  ddummy[0] = du;
#endif
  
  for(i = 1; i < VR; i++){
    ddummy[i] = ddummy[i-1]+4;
  }

  return(0);
}

void free_moment_field() {

  free(mo);
  free(df);
  free(du);
}
