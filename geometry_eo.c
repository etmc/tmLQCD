/* $Id$ */
/*******************************************************************************
 *
 *
 * Subroutines related to the lattice geometry
 *
 * The externally accessible function is
 *
 *   void geometry_eo(void)
 *     Computes the index arrays g_ipt, g_iup, g_idn, g_lexic2eo and g_eo2lexic
 *
 * original Version by
 * Author: Martin Luescher <luscher@mail.desy.ch>
 * Date: 24.10.2000
 *
 * Totally abused by M. Hasenbusch, now used for even-odd
 * decomposition of the lattice
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3adj.h"
/*#include "io.h"*/
#include "global.h"

void set_ipt() {
  int x0,x1,x2,x3;
#ifndef _NEW_GEOMETRY
  int ix=0;
  for (x0 = 0; x0 < T; x0++) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  g_ipt[x0][x1][x2][x3] = ix;
	  ix++;
	}
      }
    }
  }
#if ((defined PARALLELT) || (defined PARALLELXT))
  /* the time boundary */
  for (x0 = T; x0 < T+2; x0++) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  g_ipt[x0][x1][x2][x3] = ix;
	  ix++;
	}
      }
    }
  }
#endif
#if (defined PARALLELXT)
  /* the x boundary */
  for (x1 = LX; x1 < LX+2; x1++) {
    for (x0 = 0; x0 < T; x0++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  g_ipt[x0][x1][x2][x3] = ix;
	  ix++;
	}
      }
    }
  }
  /* The edges */
  for (x0 = T; x0 < T+2; x0++) {
    for (x1 = LX; x1 < LX+2; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  g_ipt[x0][x1][x2][x3] = ix;
	  ix++;
	}
      }
    }
  }
#endif

#else
  int i_even = 0,i_odd = VOLUME/2;
  int i_even2 = 0, i_odd2 = 0;
  for (x0 = 0; x0 < T; x0++) {
    for (x1 = 0; x1 < LX; x1++) {
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    g_ipt[x0][x1][x2][x3] = i_even;
	    i_even++;
	  } 
	  else {
	    g_ipt[x0][x1][x2][x3] = i_odd;
	    i_odd++;
	  }
	}
      }
    }
  }
#if ((defined PARALLELT) || (defined PARALLELXT))
  /* the time boundary */
  i_even = VOLUME; i_odd = VOLUME+(LX*LY*LZ)/2;
  i_even2 = VOLUME+(LX*LY*LZ); i_odd2 = VOLUME+3*(LX*LY*LZ)/2;
  for (x1 = 0; x1 < LX; x1++) {
    for (x2 = 0; x2 < LY; x2++) {
      for (x3 = 0; x3 < LZ; x3++) {
	if((T+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	  g_ipt[T][x1][x2][x3] = i_even;
	  i_even++;
	  g_ipt[T+1][x1][x2][x3] = i_odd2;
	  i_odd2++;
	} 
	else {
	  g_ipt[T][x1][x2][x3] = i_odd;
	  i_odd++;
	  g_ipt[T+1][x1][x2][x3] = i_even2;
	  i_even2++;
	}
      }
    }
  }

#endif
#if (defined PARALLELXT)
  /* the x boundary */
  i_even = VOLUME+2*LX*LY*LZ; i_odd = i_even + (T*LY*LZ)/2;
  i_even2 = VOLUME+2*LX*LY*LZ + (T*LY*LZ); i_odd2 = i_even2 + (T*LY*LZ)/2;
  for (x0 = 0; x0 < T; x0++) {
    for (x2 = 0; x2 < LY; x2++) {
      for (x3 = 0; x3 < LZ; x3++) {
	if((x0+LX+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	  g_ipt[x0][LX][x2][x3] = i_even;
	  i_even++;
	  g_ipt[x0][LX+1][x2][x3] = i_odd2;
	  i_odd2++;
	} 
	else {
	  g_ipt[x0][LX][x2][x3] = i_odd;
	  i_odd++;
	  g_ipt[x0][LX+1][x2][x3] = i_even2;
	  i_even2++;
	}
      }
    }
  }

  /* The edges */
  i_even = VOLUME+RAND; i_odd = i_even+(LY*LZ)/2;
  i_even2 = VOLUME+RAND + (LY*LZ); i_odd2 = i_even2 + (LY*LZ)/2;
  for (x2 = 0; x2 < LY; x2++) {
    for (x3 = 0; x3 < LZ; x3++) {
      if((T+LX+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	g_ipt[T][LX][x2][x3] = i_even;
	g_ipt[T][LX+1][x2][x3] = i_odd2;
	g_ipt[T+1][LX][x2][x3] = i_odd2 + (LY*LZ);
	g_ipt[T+1][LX+1][x2][x3] = i_even + 3*(LY*LZ);
	i_even++;
	i_odd2++;
      } 
      else {
	g_ipt[T][LX][x2][x3] = i_odd;
	g_ipt[T][LX+1][x2][x3] = i_even2;
	g_ipt[T+1][LX][x2][x3] = i_even2 + (LY*LZ);;
	g_ipt[T+1][LX+1][x2][x3] = i_odd + 3*(LY*LZ);
	i_odd++;
	i_even2++;
      }
    }
  }
#endif

#endif
}

void geometry() {
  
  int x0,x1,x2,x3,ix=0;
  int y0, y1, z0, z1;
  int i_even = 0,i_odd = VOLUME/2;
  int startvaluet = 0;
  int startvaluex = 0;
  int xeven[VOLUMEPLUSRAND] ALIGN;

#if (defined PARALLELT || defined PARALLELXT)
  startvaluet = 1;
#endif
#ifdef PARALLELXT
  startvaluex = 1;
#endif

  set_ipt();

  /* extended for neighbour slices at x0=-1 and x0=T */
  for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  y0 = x0;
	  y1 = x1;
	  if(x0 == -1) y0 = T+1;
	  if(x1 == -1) y1 = LX+1;
	  z0 = y0 - 1;
	  z1 = y1 - 1;
	  if(x0 == 0) z0 = T+1;
	  if(x1 == 0) z1 = LX+1;
	  ix = g_ipt[y0][y1][x2][x3];
	  /* g_proc_id*T is added to allow for odd T when the number of 
	     nodes is even */
	  if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    xeven[ix]=1;
	  } 
	  else {
	    xeven[ix]=0; 
	  }

#if ((defined PARALLELT) || (defined PARALLELXT))
	  if(x0 <  T) g_iup[ix][0] = g_ipt[x0+1][y1][x2][x3];
	  else g_iup[ix][0] = 0;
	  if(x0 > -1) g_idn[ix][0] = g_ipt[z0  ][y1][x2][x3];
	  else g_idn[ix][0] = 0;
#else
	  g_iup[ix][0] = g_ipt[(x0+T+1)%T][x1][x2][x3];
	  g_idn[ix][0] = g_ipt[(x0+t-1)%T][x1][x2][x3];
#endif
#if (defined PARALLELXT)
	  if(x1 < LX) g_iup[ix][1] = g_ipt[y0][x1+1][x2][x3];
	  else g_iup[ix][1] = 0;
	  if(x1 > -1) g_idn[ix][1] = g_ipt[y0][z1  ][x2][x3];
	  else g_idn[ix][1] = 0;
#else
	  g_iup[ix][1] = g_ipt[y0][(x1+LX+1)%LX][x2][x3];
	  g_idn[ix][1] = g_ipt[y0][(x1+LX-1)%LX][x2][x3];
#endif

	  g_iup[ix][2] = g_ipt[y0][y1][(x2+LY+1)%LY][x3];
	  g_idn[ix][2] = g_ipt[y0][y1][(x2+LY-1)%LY][x3];

	  g_iup[ix][3] = g_ipt[y0][y1][x2][(x3+LY+1)%LZ];
	  g_idn[ix][3] = g_ipt[y0][y1][x2][(x3+LY-1)%LZ];

	  if(g_proc_id == 0) printf("%d %d %d %d: %d %d  %d %d %d %d %d %d %d %d\n", x0, x1, x2, x3, ix, xeven[ix],g_iup[ix][0], g_idn[ix][0], g_iup[ix][1], g_idn[ix][1], g_iup[ix][2], g_idn[ix][2], g_iup[ix][3], g_idn[ix][3]);

	}
      }
    }
  }
  i_even=0;
  i_odd=0;
  for (ix = 0; ix < (VOLUME+RAND); ix++){
    if(xeven[ix]==1){
      g_lexic2eo[ix] = i_even;
      g_lexic2eosub[ix] = i_even;
      g_eo2lexic[i_even] = ix;
      i_even++;
    }
    else{
      g_lexic2eo[ix] = (VOLUME+RAND)/2+i_odd;
      g_lexic2eosub[ix] = i_odd;
      g_eo2lexic[(VOLUME+RAND)/2+i_odd] = ix;
      i_odd++;
    }
  }
}



