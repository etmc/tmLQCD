#include "dirty_shameful_business.h"

void ohnohack_remap_g_gauge_field(su3_tuple *gf)
{
  g_gauge_field[0] = (su3*)gf[0];
  for(i = 1; i < V; i++)
    g_gauge_field[i] = g_gauge_field[i - 1] + 4;
}

void ohnohack_remap_g_gauge_field_copy(su3_tuple *gf)
{
#ifdef _USE_HALFSPINOR
EXTERN su3 *** g_gauge_field_copy;
#elif (defined _USE_TSPLITPAR )
EXTERN su3 ** g_gauge_field_copyt;
EXTERN su3 ** g_gauge_field_copys;
#else
  g_gauge_field_copy[0] = (su3*)gf[0];
  for(i = 1; i < (VOLUME+RAND); ++i)
    g_gauge_field_copy[i] = g_gauge_field_copy[i - 1] + 8; /* FIXME Hmm, not 4. Check how this comes in!. */
#endif
}

/** What follows is copied from init_gauge_field.c for reference and should be removed! */

#  if defined _USE_HALFSPINOR
  if(back == 1) {
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
  if(back == 1) {
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
  if(back == 1) {
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
