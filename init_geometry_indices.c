/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "init_geometry_indices.h"

int *iup = NULL, *idn = NULL, *ipt = NULL, **ipt_ = NULL, ***ipt__ = NULL;

int init_geometry_indices(const int V) {

  g_idn= calloc(VOLUMEPLUSRAND, sizeof(int*));
  if(errno == ENOMEM) return(1);
  g_iup = calloc(VOLUMEPLUSRAND, sizeof(int*));
  if(errno == ENOMEM) return(2);
  ipt_ = calloc((T+2)*(LX+2)*LY, sizeof(int*));
  if(errno == ENOMEM) return(3);
  ipt__ = calloc ((T+2)*(LX+2), sizeof(int*));
  if(errno == ENOMEM) return(4);
  g_ipt = calloc(T,sizeof(int*));
  if(errno == ENOMEM) return(5);
  idn = calloc(4*V, sizeof(int));
  if(errno == ENOMEM) return(6);
  iup = calloc(4*V, sizeof(int));
  if(errno == ENOMEM) return(7);
  ipt = calloc((T+2)*(LX+2)*LY*LZ, sizeof(int));
  if(errno == ENOMEM) return(8);
  g_lexic2eo = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(9);
  g_lexic2eosub = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(10);
  g_eo2lexic = calloc(V, sizeof(int));
  if(errno == ENOMEM) return(11);

  g_idn[0] = idn;
  g_iup[0] = iup;
  ipt_[0] = ipt;
  
  ipt__[0] = ipt_;
  g_ipt[0] = ipt__;
  for(i = 1; i < V; i++){
    g_idn[i] = g_idn[i-1]+4;
    g_iup[i] = g_iup[i-1]+4;
  }
  for(i = 1; i < V; i++){
    g_gauge_field[i] = g_gauge_field[i-1]+4;
  }
  for(i = 1; i < (T+2)*(LX+2)*LY; i++){
    ipt_[i] = ipt_[i-1]+LZ;
  }
  for(i = 1; i < (T+2)*(LX+2); i++){
    ipt__[i] = ipt__[i-1]+LY;
  }
  for(i = 1; i < (T+2); i++){
    g_ipt[i] = g_ipt[i-1]+(LX+2);
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
}
