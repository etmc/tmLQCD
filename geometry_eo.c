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


#ifndef _NEW_GEOMETRY
int Index(const int x0, const int x1, const int x2, const int x3) {
  int y0, y1, y2, y3, ix;

  y0 = (x0 + T ) % T; 
  y1 = (x1 + LX) % LX; 
  y2 = (x2 + LY) % LY; 
  y3 = (x3 + LZ) % LZ;
  ix = ((y0*LX + y1)*LY + y2)*LZ + y3;
  
  y0=x0;
#if ((defined PARALLELT) || (defined PARALLELXT))
  if(x0 == T) {
    ix = VOLUME + y3 + LZ*y2 + LZ*LY*y1;
  }
  /* the slice at time -1 is put to T+1 */
  else if(x0 == -1) {
    ix = VOLUME + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
  }
#endif
#if (defined PARALLELXT)
  if(x1 == LX){
    ix = VOLUME + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }
  if(x1 == -1){
    ix = VOLUME + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }   
  /* The edges */
  if(x0 == T){
    if(x1 == LX){
      ix = VOLUME+RAND+y2*LZ+y3;
    }
    if(x1 == -1){
      ix = VOLUME+RAND+LY*LZ+y2*LZ+y3;
    }
  }
  if(x0 == -1){
    if(x1 == LX){
      ix = VOLUME+RAND+2*LY*LZ+y2*LZ+y3;
    }
    if(x1 == -1){
      ix = VOLUME+RAND+3*LY*LZ+y2*LZ+y3;
    }
  }
#endif
  /* The DBW2 stuff --> second boundary slice */
  /* This we put a the very end.              */
#if ((defined PARALLELT) || (defined PARALLELXT))
  if(x0 == T+1) { 
    ix = VOLUMEPLUSRAND + y3 + LZ*y2 + LZ*LY*y1;
#ifdef PARALLELXT
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 1*LY*LZ
	+ y2*LZ + y3;
    }
#endif
  }
  /* the slice at time -2 is put behind the one at time T+1 */
  else if(x0 == -2) {
    ix = VOLUMEPLUSRAND + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
#ifdef PARALLELXT
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 3*LY*LZ
	+ y2*LZ + y3;
    }
#endif
  }  
#endif
#if (defined PARALLELXT)
  if(x1 == LX+1) {
    if((x0 < T) && (x0 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 6*LY*LZ
	+ y2*LZ + y3;
    }
  }
  if(x1 == -2) {
    if((x0 < T) && (x0 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 5*LY*LZ
	+ y2*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 7*LY*LZ
	+ y2*LZ + y3;
    }
  }   
#endif
  return(ix);
}

#else

int Index(const int x0, const int x1, const int x2, const int x3)
{
   int y0, y1, y2, y3, ix, bndt=0, bndx=0, odd;

#if ((defined PARALLELT) || (defined PARALLELXT))
   y0 = x0;
   /* the slice at time -1 is put to T+1 */
   if(x0 == -1) y0 = T+1;
   if(x0 == -1 || x0 == T) bndt = 1;
   if(x0 == -2) y0 = 1;
   if(x0 == T+1) y0 = 0;
   if(x0 == -2 || x0 == T+1) bndt = 2;
#else
   y0 = (x0+T) % T;
#endif

#ifdef PARALLELXT 
   y1 = x1;
   /* the slice at x -1 is put to LX+1 */
   if(x1 == -1) y1=LX+1;
   if(x1 == -1 || x1 == LX) bndx = 1;
   /* the slice at x -1 is put to LX+2 */
   if(x1 == -2) y1 = 1;
   if(x1 == LX+1) y1 = 0;
   if(x1 == -2 || x1 == LX+1) bndx = 2;
#else
   y1 = (x1+LX) % LX;
#endif

   y2=(x2+LY)%LY; 
   
   y3=(x3+LZ)%LZ; 

   if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
     odd = 0;
   } 
   else {
     odd = 1;
   }

   if(bndt == 0 && bndx == 0) {
     ix = (y3 + LZ*y2 + LY*LZ*y1 + LX*LY*LZ*y0)/2 + (odd*(VOLUME))/2;
   }
   else if(bndt == 1 && bndx == 0) {
     ix = y0*LX*LY*LZ+(y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   else if(bndt == 0 && bndx == 1) {
     ix = VOLUME + 2*LX*LY*LZ 
       + (y1-LX)*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   else if(bndt == 1 && bndx == 1) {
     ix = VOLUME + RAND 
       + ((y0-T)*2 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 2 && bndx == 0) {
     ix = VOLUMEPLUSRAND 
       + y0*LX*LY*LZ + (y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   else if(bndt == 0 && bndx == 2) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ 
       + y1*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   else if(bndt == 2 && bndx == 1) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ 
       + (2*y0 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 1 && bndx == 2) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ 
       + ((y0-T)*2 + y1)*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 2 && bndx == 2) {
     printf("Should not happen in index routine!\n");
     printf("%d %d %d %d\n", x0, x1, x2 ,x3);
     ix = -1;
   }
   else {
     printf("Error in index routine!\n");
     exit(1);
   }
   return( ix );
}

#endif

void geometry(){
  
  int x0,x1,x2,x3,ix;
  int i_even,i_odd;
  int startvaluet = 0;
  int startvaluex = 0;
  int * xeven;
  
  xeven = malloc(VOLUMEPLUSRAND*sizeof(int));

#if (defined PARALLELT || defined PARALLELXT)
  startvaluet = 1;
#endif
#ifdef PARALLELXT
  startvaluex = 1;
#endif

  /* extended for neighbour slices at x0=-1 and x0=T */
  for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = 0; x2 < LY; x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=Index(x0, x1, x2, x3);

	  /* g_proc_id*T is added to allow for odd T when the number of 
	     nodes is even */
	  if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    xeven[ix]=1;
	  } 
	  else {
	    xeven[ix]=0; 
	  }

	  if(x0 >= 0 && x1 >=0) g_ipt[x0][x1][x2][x3] = ix;
	  else if(x0 < 0 && x1 < 0) {
	    g_ipt[T+1][LX+1][x2][x3] = ix;
	  }
	  else if(x0 < 0) {
	    g_ipt[T+1][x1][x2][x3] = ix;
	  }
	  else if(x1 < 0) {
	    g_ipt[x0][LX+1][x2][x3] = ix;
	  }

	  g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

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
#if (defined PARALLELT || defined PARALLELXT)
  if(g_dbw2rand != 0) {
    if(g_proc_id == 0) {printf("DBW2 stuff\n");fflush(stdout);fflush(stdout);}
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x0 = -2;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  
	  g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  g_idn[ix][0] = -1;

	  if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  x0 = T+1;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  g_iup[ix][0] = -1;
	  g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	}
      }
    }    
#ifdef PARALLELXT
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x2 = 0; x2 < LY; x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  x1 = -2;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  g_idn[ix][1] = -1;

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  x1 = LX+1;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  g_iup[ix][1] = -1;
	  g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	}
      }
    }
#endif
  }
#endif
}

static char const rcsid[] = "$Id$";


