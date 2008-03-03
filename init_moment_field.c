/* $Id$ */

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

su3adj * mo=NULL, *df=NULL, *du=NULL;

int init_moment_field(const int V, const int VR) {
  int i = 0;

/*   posix_memalign(void **memptr, size_t alignment, size_t size) */
  if( (int*)(mo = (su3adj*)calloc(4*V+1, sizeof(su3adj))) == NULL){ 
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(moment = (su3adj**)calloc(V,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
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

  if((void*)(df = (su3adj*)calloc(4*VR+1, sizeof(su3adj))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(3);
  }
  if((void*)(df0 = (su3adj**)calloc(VR,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(4);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  df0[0] = (su3adj*)(((unsigned long int)(df)+ALIGN_BASE)&~ALIGN_BASE);
#else
  df0[0] = df;
#endif
  
  for(i = 1; i < VR; i++) {
    df0[i] = df0[i-1]+4; 
  }
  
  if((void*)(du = (su3adj*)calloc(4*VR+1, sizeof(su3adj))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(5);
  }
  if((void*)(ddummy = (su3adj**)calloc(VR,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(6);
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
