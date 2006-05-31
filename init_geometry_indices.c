/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "init_geometry_indices.h"

int *iup = NULL, *idn = NULL, *ipt = NULL, **ipt_ = NULL, ***ipt__ = NULL;
int *halfpt = NULL;

int init_geometry_indices(const int V) {
  int i = 0;

  g_idn= calloc(V, sizeof(int*));
  if(errno == ENOMEM) return(1);
  g_iup = calloc(V, sizeof(int*));
  if(errno == ENOMEM) return(2);

  idn = calloc(4*V, sizeof(int));
  if(errno == ENOMEM) return(6);
  iup = calloc(4*V, sizeof(int));
  if(errno == ENOMEM) return(7);

  g_ipt = calloc(T+4,sizeof(int*));
  if(errno == ENOMEM) return(5);
  ipt__ = calloc ((T+4)*(LX+4), sizeof(int*));
  if(errno == ENOMEM) return(4);
  ipt_ = calloc((T+4)*(LX+4)*(LY+4), sizeof(int*));
  if(errno == ENOMEM) return(3);
  ipt = calloc((T+4)*(LX+4)*(LY+4)*(LZ+4), sizeof(int));
  if(errno == ENOMEM) return(8);

  g_lexic2eo = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(9);
  g_lexic2eosub = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(10);
  g_eo2lexic = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(11);

#if defined PARALLELXYZT
  g_field_z_ipt_even = calloc(T*LX*LY, sizeof(int));
  if(errno == ENOMEM) return(12);
  g_field_z_ipt_odd  = calloc(T*LX*LY, sizeof(int));
  if(errno == ENOMEM) return(13);
#endif

  halfpt = (int*)calloc(4*(T+2)*(LX+2)*(LY+2)*(LZ+2)/2, sizeof(int));
  if(errno == ENOMEM) return(14);
  g_halfpt = (int**)calloc((T+2)*(LX+2)*(LY+2)*(LZ+2)/2, sizeof(int*));
  if(errno == ENOMEM) return(15);

  g_idn[0] = idn;
  g_iup[0] = iup;

  ipt_[0] = ipt;
  ipt__[0] = ipt_;
  g_ipt[0] = ipt__;
  g_halfpt[0] = halfpt;
  for(i = 1; i < V; i++){
    g_idn[i] = g_idn[i-1]+4;
    g_iup[i] = g_iup[i-1]+4;
  }
  for(i = 1; i < (T+4)*(LX+4)*(LY+4); i++){
    ipt_[i] = ipt_[i-1]+(LZ+4);
  }
  for(i = 1; i < (T+4)*(LX+4); i++){
    ipt__[i] = ipt__[i-1]+(LY+4);
  }
  for(i = 1; i < (T+4); i++){
    g_ipt[i] = g_ipt[i-1]+(LX+4);
  }
  for(i = 1; i < (T+2)*(LX+2)*(LY+2)*(LZ+2)/2; i++) {
    g_halfpt[i] = g_halfpt[i-1] + 4;
  }
  return(0);
}

void free_geometry_indices() {
  free(idn); 
  free(iup); 
  free(ipt);
  free(ipt_);
  free(ipt__);
  free(g_ipt);
  free(g_idn);
  free(g_iup);
  free(g_eo2lexic);
  free(g_lexic2eosub);
  free(g_lexic2eo);
  free(g_halfpt);
  free(halfpt);
#ifdef PARALLELXYZT
  free(g_field_z_ipt_odd);
  free(g_field_z_ipt_even);
#endif
}
