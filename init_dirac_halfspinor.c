/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "su3.h"
#include "Hopping_Matrix.h"
#include "init_dirac_halfspinor.h"

halfspinor ** NBPointer_;

int init_dirac_halfspinor() {
  int ieo=0, i=0, j=0, k;
  int x, y, z, t, mu;
  
  NBPointer = (halfspinor***) malloc(4*sizeof(halfspinor**));
  NBPointer_ = (halfspinor**) malloc(16*(VOLUME+RAND)*sizeof(halfspinor*));
  NBPointer[0] = NBPointer_;
  NBPointer[1] = NBPointer_ + (8*(VOLUME+RAND)/2);
  NBPointer[2] = NBPointer_ + (16*(VOLUME+RAND)/2);
  NBPointer[3] = NBPointer_ + (24*(VOLUME+RAND)/2);

  /* Memory for HalfSpinor must be allocated already */
  if(HalfSpinor == NULL) {
    return(-1);    
  }

  for(ieo = 0; ieo < 2; ieo++) {
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ieo*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 4; mu++) {
	NBPointer[ieo][8*i + 2*mu + 0] = &HalfSpinor[ 8*g_lexic2eosub[g_idn[j][mu]] + 2*mu + 1];
	NBPointer[ieo][8*i + 2*mu + 1] = &HalfSpinor[ 8*g_lexic2eosub[g_iup[j][mu]] + 2*mu + 0];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i] = &HalfSpinor[ k ];
      }
      if(t == T-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i + 1] = &HalfSpinor[ k ];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 2] = &HalfSpinor[ k ];
      }
      if(x == LX-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 3] = &HalfSpinor[ k ];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 4] = &HalfSpinor[ k ];
      }
      if(y == T-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 5] = &HalfSpinor[ k ];
      }
#endif
#if (defined PARALLELXYZT)
      if(z == 0) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 6] = &HalfSpinor[ k ];
      }
      if(z == T-1) {
	k = 8*VOLUME/2 + (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 7] = &HalfSpinor[ k ];
      }
#endif
    }
    NBPointer[ieo][4*VOLUME] = NBPointer[ieo][0]; 
  }
  for(ieo = 2; ieo < 4; ieo++) {
    k=0;
    for(i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + (ieo-2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 4; mu++) {
	NBPointer[ieo][k + 2*mu]     = &HalfSpinor[k + 2*mu];
	NBPointer[ieo][k + 2*mu + 1] = &HalfSpinor[k + 2*mu + 1];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	NBPointer[ieo][k + 1] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_idn[j][0]] - VOLUME/2)];
      }
      if(t == T-1) {
	NBPointer[ieo][k] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_iup[j][0]] - VOLUME/2) ];
      }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == 0) {
	NBPointer[ieo][k + 3] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_idn[j][1]] - VOLUME/2)];
      }
      if(x == LX-1) { 
	NBPointer[ieo][k + 2] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_iup[j][1]] - VOLUME/2) ];
      }
#endif
#if ((defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == 0) {
	NBPointer[ieo][k + 5] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_idn[j][2]] - VOLUME/2)];
      }
      if(y == T-1) {
	NBPointer[ieo][k + 4] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_iup[j][2]] - VOLUME/2) ];
      }
#endif
#if (defined PARALLELXYZT)
      if(z == 0) {
	NBPointer[ieo][k + 7] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_idn[j][3]] - VOLUME/2)];
      }
      if(z == T-1) {
	NBPointer[ieo][k + 6] = &HalfSpinor[ 4*VOLUME + RAND/2 + (g_lexic2eosub[g_iup[j][3]] - VOLUME/2) ];
      }
#endif
      k+=8;
    }
    NBPointer[ieo][4*VOLUME] = NBPointer[ieo][0]; 
  }
  return(0);
}
