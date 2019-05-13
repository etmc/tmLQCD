/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

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
#include "expo.h"

su3 * gauge_field = NULL;
su3_32 * gauge_field_32 = NULL;
#ifdef _USE_TSPLITPAR
su3 * gauge_field_copyt = NULL;
su3 * gauge_field_copys = NULL;
#else
su3 * gauge_field_copy = NULL;
su3_32 * gauge_field_copy_32 = NULL;
#endif

int init_gauge_field(const int V, const int back) {
  int i=0;

#ifdef _USE_TSPLITPAR
  g_gauge_field_copyt = NULL;
  g_gauge_field_copys = NULL;
#else
  g_gauge_field_copy = NULL;
#endif

  if (g_exposu3_no_c == 0) init_exposu3();

  if((void*)(g_gauge_field = (su3**)calloc(V, sizeof(su3*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(gauge_field = (su3*)calloc(4*V+1, sizeof(su3))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  g_gauge_field[0] = (su3*)(((unsigned long int)(gauge_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_gauge_field[0] = gauge_field;
#endif
  for(i = 1; i < V; i++){
    g_gauge_field[i] = g_gauge_field[i-1]+4;
  }

#  if defined _USE_HALFSPINOR
  if(back == 1 && !lowmem_flag) {
    /*
      g_gauge_field_copy[ieo][PM][sites/2][mu]
    */
    if((void*)(g_gauge_field_copy = (su3***)calloc(2, sizeof(su3**))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(g_gauge_field_copy[0] = (su3**)calloc(VOLUME, sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    g_gauge_field_copy[1] = g_gauge_field_copy[0] + (VOLUME)/2;
    if((void*)(gauge_field_copy = (su3*)calloc(4*(VOLUME)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
    }
#    if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field_copy[0][0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#    else
    g_gauge_field_copy[0][0] = gauge_field_copy;
#    endif
    for(i = 1; i < (VOLUME)/2; i++) {
      g_gauge_field_copy[0][i] = g_gauge_field_copy[0][i-1]+4;
    }
    g_gauge_field_copy[1][0] = g_gauge_field_copy[0][0] + 2*VOLUME; 
    for(i = 1; i < (VOLUME)/2; i++) {
      g_gauge_field_copy[1][i] = g_gauge_field_copy[1][i-1]+4;
    }
  }
#  elif defined _USE_TSPLITPAR
  if(back == 1 && !lowmem_flag) {
    if((void*)(g_gauge_field_copyt = (su3**)calloc((VOLUME+RAND), sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(g_gauge_field_copys = (su3**)calloc((VOLUME+RAND), sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(gauge_field_copyt = (su3*)calloc(2*(VOLUME+RAND)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
    }
    if((void*)(gauge_field_copys = (su3*)calloc(6*(VOLUME+RAND)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
    }
#    if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field_copyt[0] = (su3*)(((unsigned long int)(gauge_field_copyt)+ALIGN_BASE)&~ALIGN_BASE);
    g_gauge_field_copys[0] = (su3*)(((unsigned long int)(gauge_field_copys)+ALIGN_BASE)&~ALIGN_BASE);
#    else
    g_gauge_field_copyt[0] = gauge_field_copyt;
    g_gauge_field_copys[0] = gauge_field_copys;
#    endif
    for(i = 1; i < (VOLUME+RAND); i++) {
      g_gauge_field_copyt[i] = g_gauge_field_copyt[i-1]+2;
      g_gauge_field_copys[i] = g_gauge_field_copys[i-1]+6;
    }
  }
#  else  /* than _USE_HALFSPINOR or _USE_TSPLITPAR */
  if(back == 1 && !lowmem_flag) {
    if((void*)(g_gauge_field_copy = (su3**)calloc((VOLUME+RAND), sizeof(su3*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(3);
    }
    if((void*)(gauge_field_copy = (su3*)calloc(8*(VOLUME+RAND)+1, sizeof(su3))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(4);
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
#  endif
  g_update_gauge_copy = 1;
  return(0);
}

void free_gauge_field() {
  free(gauge_field);
  free(g_gauge_field);
  if(!lowmem_flag){
#  if defined _USE_TSPLITPAR
    free(gauge_field_copys);
    free(gauge_field_copyt);
#  else
    free(gauge_field_copy);
#  endif
  }
}



int init_gauge_field_32(const int V, const int back) {
  if(!lowmem_flag){
    int i=0;

    g_gauge_field_copy_32 = NULL;


    if((void*)(g_gauge_field_32 = (su3_32**)calloc(V, sizeof(su3_32*))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(1);
    }
    if((void*)(gauge_field_32 = (su3_32*)calloc(4*V+1, sizeof(su3_32))) == NULL) {
      printf ("malloc errno : %d\n",errno); 
      errno = 0;
      return(2);
    }

    /*doing alignment no matter what*/
    g_gauge_field_32[0] = (su3_32*)(((unsigned long int)(gauge_field_32)+ALIGN_BASE32)&~ALIGN_BASE32);

    for(i = 1; i < V; i++){
      g_gauge_field_32[i] = g_gauge_field_32[i-1]+4;
    }

#    if defined _USE_HALFSPINOR
    if(back == 1) {
      /*
        g_gauge_field_copy[ieo][PM][sites/2][mu]
      */
      if((void*)(g_gauge_field_copy_32 = (su3_32***)calloc(2, sizeof(su3_32**))) == NULL) {
        printf ("malloc errno : %d\n",errno); 
        errno = 0;
        return(3);
      }
      if((void*)(g_gauge_field_copy_32[0] = (su3_32**)calloc(VOLUME, sizeof(su3_32*))) == NULL) {
        printf ("malloc errno : %d\n",errno); 
        errno = 0;
        return(3);
      }
      g_gauge_field_copy_32[1] = g_gauge_field_copy_32[0] + (VOLUME)/2;
      if((void*)(gauge_field_copy_32 = (su3_32*)calloc(4*(VOLUME)+1, sizeof(su3_32))) == NULL) {
        printf ("malloc errno : %d\n",errno); 
        errno = 0;
        return(4);
      }
      /* doing alignment no matter what */
      g_gauge_field_copy_32[0][0] = (su3_32*)(((unsigned long int)(gauge_field_copy_32)+ALIGN_BASE32)&~ALIGN_BASE32);

      for(i = 1; i < (VOLUME)/2; i++) {
        g_gauge_field_copy_32[0][i] = g_gauge_field_copy_32[0][i-1]+4;
      }
      g_gauge_field_copy_32[1][0] = g_gauge_field_copy_32[0][0] + 2*VOLUME; 
      for(i = 1; i < (VOLUME)/2; i++) {
        g_gauge_field_copy_32[1][i] = g_gauge_field_copy_32[1][i-1]+4;
      }
    }
#    else  /* than _USE_HALFSPINOR  */
    if(back == 1) {
      if((void*)(g_gauge_field_copy_32 = (su3_32**)calloc((VOLUME+RAND), sizeof(su3_32*))) == NULL) {
        printf ("malloc errno : %d\n",errno); 
        errno = 0;
        return(3);
      }
      if((void*)(gauge_field_copy_32 = (su3_32*)calloc(8*(VOLUME+RAND)+1, sizeof(su3_32))) == NULL) {
        printf ("malloc errno : %d\n",errno); 
        errno = 0;
        return(4);
      }

      /* doing alignment no matter what */
      g_gauge_field_copy_32[0] = (su3_32*)(((unsigned long int)(gauge_field_copy_32)+ALIGN_BASE32)&~ALIGN_BASE32);

      for(i = 1; i < (VOLUME+RAND); i++) {
        g_gauge_field_copy_32[i] = g_gauge_field_copy_32[i-1]+8;
      }
    }
#    endif
    g_update_gauge_copy_32 = 1;
  }
  return(0);
}

void free_gauge_field_32() {
  if(!lowmem_flag){
    free(gauge_field_32);
    free(g_gauge_field_32);
    free(gauge_field_copy_32);
  }
}


void convert_32_gauge_field( su3_32** gf32, su3** gf, int V){
 int i,mu;   
  for(i = 0; i < V; i++) {
    for(mu =0; mu<4; mu++){
     gf32[i][mu].c00 = (_Complex float) gf[i][mu].c00;
     gf32[i][mu].c01 = (_Complex float) gf[i][mu].c01;
     gf32[i][mu].c02 = (_Complex float) gf[i][mu].c02;
     
     gf32[i][mu].c10 = (_Complex float) gf[i][mu].c10;
     gf32[i][mu].c11 = (_Complex float) gf[i][mu].c11;
     gf32[i][mu].c12 = (_Complex float) gf[i][mu].c12;    

     gf32[i][mu].c20 = (_Complex float) gf[i][mu].c20;
     gf32[i][mu].c21 = (_Complex float) gf[i][mu].c21;
     gf32[i][mu].c22 = (_Complex float) gf[i][mu].c22;        
    }
  }
#if defined _USE_HALFSPINOR
  
  
  
  
#endif
  
}






