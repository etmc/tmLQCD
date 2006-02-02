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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"


#ifndef _NEW_GEOMETRY
int Index(const int x0, const int x1, const int x2, const int x3) {
  int y0, y1, y2, y3, ix;

  y0 = (x0 + T ) % T; 
  y1 = (x1 + LX) % LX; 
  y2 = (x2 + LY) % LY; 
  y3 = (x3 + LZ) % LZ;
  ix = ((y0*LX + y1)*LY + y2)*LZ + y3;
  
  y0=x0;
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT))
  if(x0 == T) {
    ix = VOLUME + y3 + LZ*y2 + LZ*LY*y1;
  }
  /* the slice at time -1 is put to T+1 */
  else if(x0 == -1) {
    ix = VOLUME + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
  }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT))
  if(x1 == LX){
    ix = VOLUME + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }
  if(x1 == -1){
    ix = VOLUME + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }   
  /* The edges */
  /* xt-edge */
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
  /* endif of PARALLELXT || PARALLELXYT*/
#endif

#if (defined PARALLELXYT)
  /* y-Rand */
  if(x2 == LY) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
  }
  if(x2 == -1) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;
  }
  /* the edges */
  /* yx-edge  */
  if(x1 == LX) {
    if(x2 == LY) {
      ix = VOLUME + RAND +  4*LY*LZ + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + T*LZ + y0*LZ + y3;
    }
  }
  if(x1 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND +  4*LY*LZ + 2*T*LZ + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 3*T*LZ + y0*LZ + y3;
    }
  }
  /* ty-edge */
  /* Be carefully here! Here we need y first, then t */
  /* this is because the chain is first t dir, then y direction */
  /* this is oposit to the other edges ! */
  if(x2 == LY) {
    if(x0 == T) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + y1*LZ + y3;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + LX*LZ + y1*LZ + y3;
    }
  }
  if(x2 == -1) {
    if(x0 == T) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + 2*LX*LZ + y1*LZ + y3;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND +  4*LY*LZ + 4*T*LZ + 3*LX*LZ + y1*LZ + y3;
    }
  }

  /* endif of PARALLELXYT */
#endif

  /* The DBW2 stuff --> second boundary slice */
  /* This we put a the very end.              */
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT))
  if(x0 == T+1) { 
    ix = VOLUMEPLUSRAND + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 1*LY*LZ
	+ y2*LZ + y3;
    }
# endif
# if defined PARALLELXYT
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ
	+ y1*LZ + y3;
    }    
# endif
  }
  /* the slice at time -2 is put behind the one at time T+1 */
  else if(x0 == -2) {
    ix = VOLUMEPLUSRAND + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 2*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 3*LY*LZ
	+ y2*LZ + y3;
    }
# endif
# if defined PARALLELXYT
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + LZ*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 3*LX*LZ
	+ y1*LZ + y3;
    }    
# endif
  }  
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT))
  if(x1 == LX+1) {
    if((x0 < T) && (x0 > -1) && (x2 < LY) && (x2 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    /* x2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 4*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 6*LY*LZ
	+ y2*LZ + y3;
    }
# ifdef PARALLELXYT
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 3*T*LZ 
	+ y0*LZ + y3;
    }
# endif
  }
  if(x1 == -2) {
    if((x0 < T) && (x0 > -1) && (x2 < LY) && (x2 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    }
    /* x2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 5*LY*LZ
	+ y2*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 7*LY*LZ
	+ y2*LZ + y3;
    }
# ifdef PARALLELXYT
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 1*T*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 3*T*LZ
	+ y0*LZ + y3;
    }
# endif
  }   
#endif
#if defined PARALLELXYT
  if(x2 == LY+1) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
    }
    /* y2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ
	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ
	+ y0*LZ + y3;
    }
    /* y2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 5*LX*LZ
	+ y2*LZ + y3;
    }
  }
  if(x2 == -2) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;      
    }
    /* y2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 5*T*LZ
	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 7*T*LZ
	+ y0*LZ + y3;
    }
    /* y2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 7*LX*LZ
	+ y2*LZ + y3;
    }
  }
#endif
  return(ix);
}

#else
/* now for even/odd geometry in the gauge fields. */

int Index(const int x0, const int x1, const int x2, const int x3)
{
  int y0, y1, y2, y3, ix, bndt=0, bndx=0, odd, bndy=0;

#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT))
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
  
#if ((defined PARALLELXT) || (defined PARALLELXYT))
  y1 = x1;
  /* the slice at x -1 is put to LX+1 */
  if(x1 == -1) y1=LX+1;
  if(x1 == -1 || x1 == LX) bndx = 1;
  /* the slice at x -2 is put to LX+2 */
  if(x1 == -2) y1 = 1;
  if(x1 == LX+1) y1 = 0;
  if(x1 == -2 || x1 == LX+1) bndx = 2;
#else
  y1 = (x1+LX) % LX;
#endif

#if defined PARALLELXYT
  y2 = x2;
  /* the slice at y -1 is put to LY+1 */
  if(x2 == -1) y2=LY+1;
  if(x2 == -1 || x2 == LY) bndy = 1;
  /* the slice at y -2 is put to LY+2 */
  if(x2 == -2) y2 = 1;
  if(x2 == LY+1) y2 = 0;
  if(x2 == -2 || x2 == LY+1) bndy = 2;
#else
   y2 = (x2+LY) % LY; 
#endif

   y3 = (x3+LZ) % LZ; 
   /* even or odd point? */
   if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX+g_proc_coords[2]*LY)%2 == 0) {
     odd = 0;
   } 
   else {
     odd = 1;
   }

   /* The local volume */
   if(bndt == 0 && bndx == 0 && bndy == 0) {
     ix = (y3 + LZ*y2 + LY*LZ*y1 + LX*LY*LZ*y0)/2 + (odd*(VOLUME))/2;
   }
   /* The time boundary */
   else if(bndt == 1 && bndx == 0 && bndy == 0) {
     ix = y0*LX*LY*LZ+(y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   /* The x boundary */
   else if(bndt == 0 && bndx == 1 && bndy == 0) {
     ix = VOLUME + 2*LX*LY*LZ 
       + (y1-LX)*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   /* the y boundary */
   else if(bndt == 0 && bndx == 0 && bndy == 1) {
     ix = VOLUME + 2*LZ*(LX*LY+T*LY) 
       + (y2-LY)*T*LX*LZ + (y0*LX*LZ + y3 + LZ*y1)/2 + (odd*(LX*LZ*T))/2;
   }
   /* the xt edges */
   else if(bndt == 1 && bndx == 1 && bndy == 0) {
     ix = VOLUME + RAND 
       + ((y0-T)*2 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   /* the xy edges */
   else if(bndt == 0 && bndx == 1 && bndy == 1) {
     ix = VOLUME + RAND + 4*LY*LZ
       + ((y1-LX)*2 + (y2-LY))*(T*LZ) + (y3 + LZ*y0)/2 + (odd*(T*LZ))/2;
   }
   /* the ty edges */
   else if(bndt == 1 && bndx == 0 && bndy == 1) {
     ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ
       + ((y0-T)*2 + (y2-LY))*(LX*LZ) + (y3 + LZ*y1)/2 + (odd*(LX*LZ))/2;
   }
   /* now comes rectangular gauge action stuff */
   else if(bndt == 2 && bndx == 0 && bndy == 0) {
     ix = VOLUMEPLUSRAND 
       + y0*LX*LY*LZ + (y3 + LZ*y2 + LY*LZ*y1)/2 + (odd*(LX*LY*LZ))/2;
   }
   else if(bndt == 0 && bndx == 2 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ 
       + y1*T*LY*LZ + (y0*LY*LZ + y3 + LZ*y2)/2 + (odd*(LY*LZ*T))/2;
   }
   else if(bndt == 2 && bndx == 1 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ 
       + (2*y0 + (y1-LX))*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 1 && bndx == 2 && bndy == 0) {
     ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 4*LY*LZ 
       + ((y0-T)*2 + y1)*(LY*LZ) + (y3 + LZ*y2)/2 + (odd*(LY*LZ))/2;
   }
   else if(bndt == 2 && bndx == 2) {
     printf("Should not happen in index routine!\n");
     printf("%d %d %d %d\n", x0, x1, x2 ,x3);
     ix = -1;
   }
   else if(bndt == 1 && bndx == 1 && bndy == 1) {
     printf("Should not be on three boundaries in index routine!\n");
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
  int y0, y1, y2, y3;
  int i_even,i_odd;
  int startvaluet = 0;
  int startvaluex = 0;
  int startvaluey = 0;
  int * xeven;
  
  xeven = malloc(VOLUMEPLUSRAND*sizeof(int));

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT)
  startvaluet = 1;
#endif
#if (defined PARALLELXT || defined PARALLELXYT)
  startvaluex = 1;
#endif
#if (defined PARALLELXYT)
  startvaluey = 1;
#endif

  /* extended for boundary slices */
  for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++){
	for (x3 = 0; x3 < LZ; x3++){
	  ix=Index(x0, x1, x2, x3);
	  y0=x0; y1=x1; y2=x2; y3=x3;
	  if(x0 == -1) {
	    y0 = T+1;
	  }
	  if(x1 == -1) {
	    y1 = LX+1;
	  }
	  if(x2 == -1) {
	    y2 = LY+1;
	  }




	  if((x0 == T  && x1 == LX && x2 == LY) ||
	     (x0 == -1 && x1 == -1 && x2 == -1) ||
	     (x0 == -1 && x1 == LX && x2 == -1) ||
	     (x0 == -1 && x1 == -1 && x2 == LY) ||
	     (x0 == T  && x1 == -1 && x2 == -1) ||
	     (x0 == -1 && x1 == LX && x2 == LY) ||
	     (x0 == T  && x1 == -1 && x2 == LY) ||
	     (x0 == T  && x1 == LX && x2 == -1)) {
	    /* Should not be needed, set it to -1 */
	    g_ipt[y0][y1][y2][y3] = -1;
	  }
	  else {
	    g_ipt[y0][y1][y2][y3] = ix;
	    /* g_proc_id*T|LX|LY is added to allow for odd T|LX|LY when the number of 
	       nodes is even */	
	    if((x0+x1+x2+x3+g_proc_coords[0]*T+g_proc_coords[1]*LX+g_proc_coords[2]*LY)%2==0) {
	      xeven[ix]=1;
	    } 
	    else {
	      xeven[ix]=0; 
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
  }
  i_even=0;
  i_odd=0;
  /*For the spinor fields we need only till VOLUME+RAND */
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

  /* The rectangular gauge action part */
  /* Everything is stored behind VOLUMEPLUSRAND-1 !*/
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT)
  if(g_dbw2rand != 0) {
    if(g_proc_id == 0) {
      printf("# Initialising rectangular gauge action stuff\n");
      fflush(stdout);
    }
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  /* t2 Rand and t2x and t2y */
	  x0 = -2;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  
	  g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  g_idn[ix][0] = -1;

	  if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

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

	  if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	}
      }
    }    
#if (defined PARALLELXT || PARALLELXYT)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  /* x2-Rand and x2t and x2y */
	  x1 = -2;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  g_idn[ix][1] = -1;

	  if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

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

	  if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	}
      }
    }
#endif
#ifdef PARALLELXYT
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++) {
	for (x3 = 0; x3 < LZ; x3++) {
	  /* y2-Rand y2t and y2x */
	  x2 = -2;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	  g_idn[ix][2] = -1;

	  g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	  g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  x2 = LY+1;
	  ix = Index(x0, x1, x2, x3);
	  if(ix < VOLUMEPLUSRAND) {
	    printf("#### %d %d %d %d\n",x0, x1, x2, x3);
	  }
	  if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	  if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	  if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	  if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	  g_iup[ix][2] = -1;
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


