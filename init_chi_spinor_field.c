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

spinor * sp_up = NULL;
spinor * sp_dn = NULL;

int init_chi_up_spinor_field(const int V, const int nr) {
  int i = 0;

  if((void*)(sp_up = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(g_chi_up_spinor_field = malloc(nr*sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_up_spinor_field[0] = (spinor*)(((unsigned long int)(sp_up)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_up_spinor_field[0] = sp_up;
#endif
  
  for(i = 1; i < nr; i++){
    g_chi_up_spinor_field[i] = g_chi_up_spinor_field[i-1]+V;
  }

  return(0);
}

void free_chi_up_spinor_field() {

  free(sp_up);
}


int init_chi_dn_spinor_field(const int V, const int nr) {
  int i = 0;

  if((void*)(sp_dn = (spinor*)calloc(nr*V+1, sizeof(spinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(g_chi_dn_spinor_field = malloc(nr*sizeof(spinor*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  g_chi_dn_spinor_field[0] = (spinor*)(((unsigned long int)(sp_dn)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_chi_dn_spinor_field[0] = sp_dn;
#endif
  
  for(i = 1; i < nr; i++){
    g_chi_dn_spinor_field[i] = g_chi_dn_spinor_field[i-1]+V;
  }

  return(0);
}

void free_chi_dn_spinor_field() {

  free(sp_dn);
}

