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

spinor * sp_cup = NULL;
spinor * sp_cdn = NULL;

int init_chi_up_copy(const int V) {
  int i = 0;

  sp_cup = (spinor*)calloc(V+1, sizeof(spinor));
  if(errno == ENOMEM) {
    return(1);
  }
  g_chi_up_copy = malloc(sizeof(spinor*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_up_copy = (spinor*)(((unsigned long int)(sp_cup)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_up_copy = sp_cup;
#endif
  
  return(0);
}

void free_chi_up_copy() {

  free(sp_cup);
}


int init_chi_dn_copy(const int V) {
  int i = 0;

  sp_cdn = (spinor*)calloc(V+1, sizeof(spinor));
  if(errno == ENOMEM) {
    return(1);
  }
  g_chi_dn_copy = malloc(sizeof(spinor*));
  if(errno == ENOMEM) {
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_dn_copy = (spinor*)(((unsigned long int)(sp_cdn)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_dn_copy = sp_cdn;
#endif

  return(0);
}

void free_chi_dn_copy() {

  free(sp_cdn);
}



