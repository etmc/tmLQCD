/* $Id$ */
/*******************************************************************************
 *
 *
 * Subroutines related to the lattice geometry
 *
 * The externally accessible function is
 *
 *   void geometry_eo(void)
 *     Computes the index arrays g_ipt, g_iup, g_idn, trans1 and trans2
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

#define OlD 

#ifdef OlD
int Index(const int x0, const int x1, const int x2, const int x3)
{
   int y0, y1, y2, y3;

   y0=x0;
#ifdef MPI
   /* the slice at time -1 is put to T+1 */
   if(x0 == -1) y0=T+1;
#else
   /* Without parallelisation */
   if(x0 == -1) y0 = T-1;
   else if(x0 == T) y0 = 0;
#endif
   y1=x1;
   if(x1==L) y1=0;
   else if(x1==-1) y1=L-1;
   y2=x2;
   if(x2==L) y2=0;
   else if(x2==-1) y2=L-1;
   y3=x3;
   if(x3==L) y3=0;
   else if(x3==-1) y3=L-1;

   return(y3 + L*y2 + L*L*y1 + L*L*L*y0);
}

#else

int Index(const int x0, const int x1, const int x2, const int x3)
{
   int y0, y1, y2, y3, ix, bnd=0;

   y0=x0;
   /* the slice at time -1 is put to T+1 */
   if(x0 == -1) y0=T+1;
   if(x0 == -1 || x0 == T) bnd = 1;

   y1=x1;
   if(x1==L) y1=0;
   else if(x1==-1) y1=L-1;

   y2=x2;
   if(x2==L) y2=0;
   else if(x2==-1) y2=L-1;

   y3=x3;
   if(x3==L) y3=0;
   else if(x3==-1) y3=L-1;
   
   if((x0+x1+x2+x3+g_proc_id*T)%2==0) {
     if(bnd ==0) {
       ix = (y3 + L*y2 + L*L*y1 + L*L*L*y0)/2;
     }
     else {
       ix = y0*L*L*L+(y3 + L*y2 + L*L*y1)/2;
     }
   }
   else {
     if(bnd == 0) {
       ix = (y3 + L*y2 + L*L*y1 + L*L*L*y0)/2+(VOLUME)/2;
     }
     else {
       ix = y0*L*L*L+(y3 + L*y2 + L*L*y1)/2+L*L*L/2;
     }
   }

   return( ix );
}

#endif

void geometry(){
  
  int x0,x1,x2,x3,ix;
  int i_even,i_odd;
  int startvalue = 1;
#ifdef MPI
  startvalue = 0;
#endif
  
  /* extended for neighbour slices at x0=-1 and x0=T */
  for (x0=-startvalue;x0<(T+startvalue);x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  ix=Index(x0,x1,x2,x3);

	  /* g_proc_id*T is added to allow for odd T when the number of 
	     nodes is even */
	  if((x0+x1+x2+x3+g_proc_id*T)%2==0){
	    xeven[ix]=1;
	  } 
	  else xeven[ix]=0; 
	  if(x0 >= 0) g_ipt[x0][x1][x2][x3] = ix;

	  g_iup[ix][0]=Index(x0+1,x1,x2,x3);
	  if(x0 >= 0) g_idn[ix][0] = Index(x0-1,x1,x2,x3);

	  g_iup[ix][1]=Index(x0,x1+1,x2,x3);
	  g_idn[ix][1]=Index(x0,x1-1,x2,x3);

	  g_iup[ix][2]=Index(x0,x1,x2+1,x3);
	  g_idn[ix][2]=Index(x0,x1,x2-1,x3);

	  g_iup[ix][3]=Index(x0,x1,x2,x3+1);
	  g_idn[ix][3]=Index(x0,x1,x2,x3-1);
	  
	}
      }
    }
  }
  i_even=0;
  i_odd=0;
  for (ix=0;ix<(VOLUME+RAND);ix++){
    if(xeven[ix]==1){
      trans1[ix]=i_even;
      trans2[i_even]=ix;
      i_even++;
    }
    else{
      trans1[ix]=(VOLUME+RAND)/2+i_odd;
      trans2[(VOLUME+RAND)/2+i_odd]=ix;
      i_odd++;
    }
  }
}



