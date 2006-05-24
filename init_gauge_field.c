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
#include "init_gauge_field.h"

su3 * gauge_field = NULL;
su3 * gauge_field_copy = NULL;

int init_gauge_field(const int V, const int back) {
  int i=0;

  g_gauge_field_copy = NULL;
  g_gauge_field = calloc(V, sizeof(su3*));
  if(errno == ENOMEM) {
    return(1);
  }
  gauge_field = calloc(4*V+1, sizeof(su3));
  if(errno == ENOMEM) {
    return(1);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  g_gauge_field[0] = (su3*)(((unsigned long int)(gauge_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_gauge_field[0] = gauge_field;
#endif
  for(i = 1; i < V; i++){
    g_gauge_field[i] = g_gauge_field[i-1]+4;
  }

  if(back == 1) {
    g_gauge_field_copy = calloc((VOLUME+RAND), sizeof(su3*));
    if(errno == ENOMEM) {
      return(2);
    }
    gauge_field_copy = calloc(8*(VOLUME+RAND)+1, sizeof(su3));
    if(errno == ENOMEM) {
      return(2);
    }
#  if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field_copy[0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#  else
    g_gauge_field_copy[0] = gauge_field_copy;
#  endif
    for(i = 1; i < (VOLUME+RAND); i++) {
      g_gauge_field_copy[i] = g_gauge_field_copy[i-1]+8;
    }
  }
  return(0);
}

void free_gauge_field() {
  free(gauge_field);
  free(g_gauge_field);
  free(gauge_field_copy);
  free(g_gauge_field_copy);
}
