/*****************************************************************************
 * Copyright (C) 2001 Martin Hasenbusch
 *               2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Modified by Jenifer Gonzalez Lopez 31.03.2009
 *
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
 *
 *
 * Subroutines related to the lattice geometry
 *
 * The externally accessible function is
 *
 *   void geometry_eo(void)
 *     Computes the index arrays g_ipt, g_iup, g_idn, g_lexic2eo and g_eo2lexic
 *
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "mpi_init.h"

void Hopping_Matrix_Indices(void);

#if ((defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ))

/* This is the version of the function Index  introduced for Aurora-like parallelizations (mainly xyz)  */
int Index(const int x0, const int x1, const int x2, const int x3) {
  /* defined for all points in the internal lattice */
  /* and for those points in the external lattice such that: */
  /* - up to 2 directions out of the lattice, */
  /* - one direction up to distance 2 out of the lattice, */ 
  /* - the other direction up to distance 1 out of the lattice. */

  int y0, y1, y2, y3, ix;

  y0 = (x0 + T ) % T; 
  y1 = (x1 + LX) % LX; 
  y2 = (x2 + LY) % LY; 
  y3 = (x3 + LZ) % LZ;
  ix = ((y0*LX + y1)*LY + y2)*LZ + y3;

  /* x-Rand */
  if(x1 == LX){
    ix = VOLUME + y0*LY*LZ + y2*LZ + y3;
  }
  if(x1 == -1){
    ix = VOLUME + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
  }   

#if (defined PARALLELXY || defined PARALLELXYZ)
  /* y-Rand */
  if(x2 == LY) {
    ix = VOLUME + 2*T*LY*LZ    + y0*LX*LZ + y1*LZ + y3;
  }
  if(x2 == -1) {
    ix = VOLUME + 2*T*LY*LZ + T*LX*LZ    + y0*LX*LZ + y1*LZ + y3;
  }
  /* yx-edge  */
  if(x1 == LX) {
    if(x2 == LY) {
      ix = VOLUME + RAND    + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + T*LZ    + y0*LZ + y3;
    }
  }
  if(x1 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 2*T*LZ    + y0*LZ + y3;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 3*T*LZ    + y0*LZ + y3;
    }
  }
#endif /* endif of PARALLELXY  || PARALLELXYZ */

#if defined PARALLELXYZ
  /* z-Rand */
  if(x3 == LZ) {
    ix = VOLUME + 2*T*LY*LZ + 2*T*LX*LZ     + y0*LX*LY + y1*LY + y2;
  }
  if(x3 == -1) {
    ix = VOLUME + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY     + y0*LX*LY + y1*LY + y2;
  }
  /* zx-edge  */
  if(x1 == LX) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*T*LZ      + y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*T*LZ + T*LY 	+ y0*LY + y2;
    }
  }
  if(x1 == -1) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*T*LZ + 2*T*LY	+ y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*T*LZ + 3*T*LY 	+ y0*LY + y2;
    }
  }
  /* zy-edge */
  if(x3 == LZ) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*T*LZ + 4*T*LY 	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*T*LZ + 4*T*LY + 2*T*LX 	+ y0*LX + y1;
    }
  }
  if(x3 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*T*LZ + 4*T*LY + T*LX 	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*T*LZ + 4*T*LY + 3*T*LX 	+ y0*LX + y1;
    }
  }

#endif /* endif of PARALLELXYZ */

  /* The DBW2 stuff --> second boundary slice */
  /* This we put a the very end.              */

  /* x2-rand+ */
  if(x1 == LX+1) {
    ix = VOLUMEPLUSRAND + y0*LY*LZ + y2*LZ + y3;
# if (defined PARALLELXY || defined PARALLELXYZ)
    /* x2y */
    if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND 	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 1*T*LZ 	+ y0*LZ + y3;
    }
# endif /* endif of PARALLELXY || PARALLELXYZ  */
# if defined PARALLELXYZ
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 4*T*LY 	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 5*T*LY 	+ y0*LY + y2;      
    }
# endif /* endif of PARALLELXYZ  */
  }
  /* x2-rand- */
  if(x1 == -2) {
    ix = VOLUMEPLUSRAND + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
# if (defined PARALLELXY || defined PARALLELXYZ)
    /* x2y */
    if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 2*T*LZ 	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 3*T*LZ	+ y0*LZ + y3;
    }
# endif /* endif of PARALLELXY || PARALLELXYZ  */
# if defined PARALLELXYZ
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 6*T*LY	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 7*T*LY	+ y0*LY + y2;      
    }
# endif /* endif of  PARALLELXYZ  */
  }   
#if (defined PARALLELXY || defined PARALLELXYZ)
  /* y2-rand+ */
  if(x2 == LY+1) {
    ix = VOLUMEPLUSRAND + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
    /* y2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 4*T*LZ	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 6*T*LZ	+ y0*LZ + y3;
    }
#  if defined PARALLELXYZ
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 4*T*LX	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 5*T*LX	+ y0*LX + y1;      
    }
#  endif /* endif of PARALLELXYZ  */
  }
  /* y2-rand- */
  if(x2 == -2) {
    ix = VOLUMEPLUSRAND + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;      
    /* y2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 5*T*LZ	+ y0*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 7*T*LZ	+ y0*LZ + y3;
    }
#  if defined PARALLELXYZ
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 6*T*LX	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 7*T*LX	+ y0*LX + y1;      
    }
# endif /* endif of PARALLELXYZ  */
  }
#endif /* endif of PARALLELXY || PARALLELXYZ  */
#if defined PARALLELXYZ
  /* z2-rand+ */
  if(x3 == LZ+1) {
    ix = VOLUMEPLUSRAND + 2*T*LY*LZ + 2*T*LX*LZ    + y0*LX*LY + y1*LY + y2;
    /* z2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ 	+ y0*LY + y2;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 2*T*LY	+ y0*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY 	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 2*T*LX	+ y0*LX + y1;      
    }
  }
  /* z2-rand- */
  if(x3 == -2) {
    ix = VOLUMEPLUSRAND + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY    + y0*LX*LY + y1*LY + y2;
    /* z2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + T*LY 	+ y0*LY + y2;
    }
    else if(x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 3*T*LY	+ y0*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 1*T*LX	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*T*LZ + 8*T*LY + 3*T*LX	+ y0*LX + y1;      
    }
  }
#endif /* endif of PARALLELXYZ  */

  return(ix);
}

#else /* original version of Index(): used for no parallelization  or PARALLEL*T */

int Index(const int x0, const int x1, const int x2, const int x3) {
  int y0, y1, y2, y3, ix;

#ifdef  WITHLAPH
  y0 = x0;
#else
  y0 = (x0 + T ) % T; 
#endif
  y1 = (x1 + LX) % LX; 
  y2 = (x2 + LY) % LY; 
  y3 = (x3 + LZ) % LZ;
  ix = ((y0*LX + y1)*LY + y2)*LZ + y3;
  
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  if(x0 == T) {
    ix = VOLUME + y3 + LZ*y2 + LZ*LY*y1;
  }
  /* the slice at time -1 is put to T+1 */
  else if(x0 == -1) {
    ix = VOLUME + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
  }
#endif
#if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
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
  
#endif /* endif of PARALLELXT || PARALLELXYT || PARALLELXYZT */

#if (defined PARALLELXYT || defined PARALLELXYZT)
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

#endif /* endif of PARALLELXYT  || PARALLELXYZT */
#if defined PARALLELXYZT
  /* z-Rand */
  if(x3 == LZ) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ 
      + y0*LX*LY + y1*LY + y2;
  }
  if(x3 == -1) {
    ix = VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY 
      + y0*LX*LY + y1*LY + y2;
  }
  /* the edges */
  /* zx-edge  */
  if(x1 == LX) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ 
	+ y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + T*LY 
	+ y0*LY + y2;
    }
  }
  if(x1 == -1) {
    if(x3 == LZ) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY 
	+ y0*LY + y2;
    }
    if(x3 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 3*T*LY 
	+ y0*LY + y2;
    }
  }
  /* tz-edge */
  if(x3 == LZ) {
    if(x0 == T) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY 
	+ y1*LY + y2;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + LX*LY 
	+ y1*LY + y2;
    }
  }
  if(x3 == -1) {
    if(x0 == T) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY 
	+ y1*LY + y2;
    }
    if(x0 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 3*LX*LY 
	+ y1*LY + y2;
    }
  }
  /* zy-edge */
  if(x3 == LZ) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY 
	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX 
	+ y0*LX + y1;
    }
  }
  if(x3 == -1) {
    if(x2 == LY) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + T*LX 
	+ y0*LX + y1;
    }
    if(x2 == -1) {
      ix = VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 3*T*LX 
	+ y0*LX + y1;
    }
  }


#endif /* endif of PARALLELXYZT */

  /* The DBW2 stuff --> second boundary slice */
  /* This we put a the very end.              */
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
  if(x0 == T+1) { 
    ix = VOLUMEPLUSRAND + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 1*LY*LZ
	+ y2*LZ + y3;
    }
# endif /* endif of PARALLELXT || PARALLELXYT || PARALLELXYZT  */
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ
	+ y1*LZ + y3;
    }    
# endif /* endif of PARALLELXYT || PARALLELXYZT  */
# if defined PARALLELXYZT
    /* t2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ
	+ y1*LY + y2;
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY
	+ y1*LY + y2;
    }
# endif /* endif of PARALLELXYZT  */
  }
  /* the slice at time -2 is put behind the one at time T+1 */
  else if(x0 == -2) {
    ix = VOLUMEPLUSRAND + LX*LY*LZ + y3 + LZ*y2 + LZ*LY*y1;
# if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    /* t2x */
    if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 2*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 3*LY*LZ
	+ y2*LZ + y3;
    }
# endif /* endif of PARALLELXT || PARALLELXYT || PARALLELXYZT  */
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* t2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + LX*LZ
	+ y1*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 3*LX*LZ
	+ y1*LZ + y3;
    }    
# endif /* endif of PARALLELXYT || PARALLELXYZT  */
# if defined PARALLELXYZT
    /* t2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + LX*LY
	+ y1*LY + y2;
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 3*LX*LY
	+ y1*LY + y2;
    }
# endif /* endif of PARALLELXYZT  */
  }  
#endif  /* endif of PARALLELT || PARALLELXT || PARALLELXYT || PARALLELXYZT  */
#if ((defined PARALLELXT) || (defined PARALLELXYT) || defined PARALLELXYZT)
  if(x1 == LX+1) {
    ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    /* x2t */
    if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 4*LY*LZ
	+ y2*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 6*LY*LZ
	+ y2*LZ + y3;
    }
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 1*T*LZ 
	+ y0*LZ + y3;
    }
# endif /* endif of PARALLELXYT || PARALLELXYZT  */
# if defined PARALLELXYZT
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY
	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 5*T*LY
	+ y0*LY + y2;      
    }
# endif /* endif of PARALLELXYZT  */
  }
  if(x1 == -2) {
    ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + T*LY*LZ + y0*LY*LZ + y2*LZ + y3;
    /* x2t */
    if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 5*LY*LZ
	+ y2*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 7*LY*LZ
	+ y2*LZ + y3;
    }
# if (defined PARALLELXYT || defined PARALLELXYZT)
    /* x2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ
	+ y0*LZ + y3;
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 3*T*LZ
	+ y0*LZ + y3;
    }
# endif /* endif of PARALLELXYT || PARALLELXYZT  */
# if defined PARALLELXYZT
    /* x2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY
	+ y0*LY + y2;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 7*T*LY
	+ y0*LY + y2;      
    }
# endif /* endif of  PARALLELXYZT  */
  }   
#endif /* endif of PARALLELXT || PARALLELXYT || PARALLELXYZT  */
#if (defined PARALLELXYT || defined PARALLELXYZT)
  if(x2 == LY+1) {
    ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + y0*LX*LZ + y1*LZ + y3;
    /* y2x */
    if(x1 == LX) {
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
	+ y1*LZ + y3;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 5*LX*LZ
	+ y1*LZ + y3;
    }
#  if defined PARALLELXYZT
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX
	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 5*T*LX
	+ y0*LX + y1;      
    }
#  endif /* endif of PARALLELXYZT  */
  }
  if(x2 == -2) {
    ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + T*LX*LZ + y0*LX*LZ + y1*LZ + y3;      
    /* y2x */
    if(x1 == LX) {
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
	+ y1*LZ + y3;
    }
    else if (x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 7*LX*LZ
	+ y1*LZ + y3;
    }
#  if defined PARALLELXYZT
    /* y2z */
    else if(x3 == LZ) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX
	+ y0*LX + y1;      
    }
    else if(x3 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 7*T*LX
	+ y0*LX + y1;      
    }
# endif /* endif of PARALLELXYZT  */
  }
#endif /* endif of PARALLELXYT || PARALLELXYZT  */
#if defined PARALLELXYZT
  /* z2-Rand */
  if(x3 == LZ+1) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x2 > -1) && (x2 < LY)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + y0*LX*LY + y1*LY + y2;
    }
    /* z2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY 
	+ y0*LY + y2;
    }
    else if (x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY
	+ y0*LY + y2;
    }
    /* z2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY
	+ y1*LY + y2;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 5*LX*LY
	+ y1*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY 
	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX
	+ y0*LX + y1;      
    }
  }
  if(x3 == -2) {
    if((x0 < T) && (x0 > -1) && (x1 < LX) && (x1 > -1) && (x2 > -1) && (x2 < LY)) {
      ix = VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY + y0*LX*LY + y1*LY + y2;
    }
    /* z2x */
    else if(x1 == LX) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + T*LY 
	+ y0*LY + y2;
    }
    else if(x1 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 3*T*LY
	+ y0*LY + y2;
    }
    /* z2t */
    else if(x0 == T) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY
	+ y1*LY + y2;
    }
    else if(x0 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 7*LX*LY
	+ y1*LY + y2;
    }
    /* z2y */
    else if(x2 == LY) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 1*T*LX
	+ y0*LX + y1;      
    }
    else if(x2 == -1) {
      ix = VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 3*T*LX
	+ y0*LX + y1;      
    }
  }
#endif /* endif of PARALLELXYZT  */
/*   if(ix == 372) { */
/*     printf("## %d %d %d %d ix = %d, %d %d %d %d\n", x0, x1, x2, x3, ix, T, LX, LY, LZ); */
/*   } */
  return(ix);
}

#endif /* PARALLEL???  */

void geometry(){
  
  int x0,x1,x2,x3,ix;
  int y0, y1, y2, y3, j;
  int bndcnt=0;
  int i_even,i_odd;
  int startvaluet = 0;
  int startvaluex = 0;
  int startvaluey = 0;
  int startvaluez = 0;
  int * xeven;
#if defined MPI
  int isp, *ones, *oneS, *oneL;
  int lsliceS, lsliceL, check_struct_zt;
#endif

  xeven = malloc(VOLUMEPLUSRAND*sizeof(int));

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  startvaluet = 1;
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  startvaluex = 1;
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  startvaluey = 1;
#endif
#if (defined PARALLELXYZT || defined PARALLELXYZ )
  startvaluez = 1;
#endif

  /* extended for boundary slices */
  for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++){
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++){
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;

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
	  if(x3 == -1) {
	    y3 = LZ+1;
	  }
	  if(bndcnt > 2) {
	    /* Should not be needed, set it to -1 */
	    g_ipt[y0][y1][y2][y3] = -1;
	  }
	  else {
	    ix=Index(x0, x1, x2, x3);
	    g_ipt[y0][y1][y2][y3] = ix;
	    /* g_proc_id*T|LX|LY|LZ is added to allow for odd T|LX|LY|LZ when the number of 
	       nodes is even */	
	    if((x0 + x1 + x2 + x3 + 
		g_proc_coords[0]*T + g_proc_coords[1]*LX + 
		g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) {
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

	    if(ix<VOLUME){
	      g_coord[ix][0]=x0+g_proc_coords[0]*T;
	      g_coord[ix][1]=x1+g_proc_coords[1]*LX;
	      g_coord[ix][2]=x2+g_proc_coords[2]*LY;
	      g_coord[ix][3]=x3+g_proc_coords[3]*LZ;
	    }
	  }
	}
      }
    }
  }


#ifdef WITHLAPH
  tempT=T;
  T=1;
  tempV=VOLUME;
  VOLUME=SPACEVOLUME;
  tempR=RAND;
  RAND=SPACERAND;
  x0=0;
  for (x1 = 0; x1 < (LX); x1++){
    for (x2 = 0; x2 < (LY); x2++){
      for (x3 = 0; x3 < (LZ); x3++){
	ix=Index(x0, x1, x2, x3);
	g_iup3d[ix][0] = -1;
	g_idn3d[ix][0] = -1;
	g_iup3d[ix][1] = Index(x0, x1+1, x2, x3);
	g_idn3d[ix][1] = Index(x0, x1-1, x2, x3);
	g_iup3d[ix][2] = Index(x0, x1, x2+1, x3);
	g_idn3d[ix][2] = Index(x0, x1, x2-1, x3);
	g_iup3d[ix][3] = Index(x0, x1, x2, x3+1);
	g_idn3d[ix][3] = Index(x0, x1, x2, x3-1);
      }
    }
  }
  T=tempT;
  VOLUME=tempV;
  RAND=tempR;
#endif

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

  for(j=0; j<4; j++){  // NEW GIUPDNEO
    for(ix = 0; ix < (VOLUME)/2; ix++){
      g_iup_eo[ix][j]=g_lexic2eosub[g_iup[g_eo2lexic[ix]][j]];
      g_idn_eo[ix][j]=g_lexic2eosub[g_idn[g_eo2lexic[ix]][j]];
    }
    for(ix = (VOLUME+RAND)/2; ix < VOLUME+RAND/2; ix++){
      g_iup_eo[ix][j]=g_lexic2eosub[g_iup[g_eo2lexic[ix]][j]];
      g_idn_eo[ix][j]=g_lexic2eosub[g_idn[g_eo2lexic[ix]][j]];
    }
  }

/* this if statement will be removed in future and _INDEX_INDEP_GEOM will be the default */
#if defined _INDEX_INDEP_GEOM 
# ifdef _USE_TSPLITPAR
  /* compute the first point (in eo system) of each timeslice */ 
  for (x0 = 0; x0 < T; x0++){
    ix = Index(x0,0,0,0);
    if(xeven[ix]==1){
      g_1st_eot[x0][0]=g_lexic2eo[ix];
      g_1st_eot[x0][1]=g_lexic2eo[ix+1];
    }else{
      g_1st_eot[x0][0]=g_lexic2eo[ix+1];
      g_1st_eot[x0][1]=g_lexic2eo[ix];
    }

    /* Starting points of ?t slices (eo) */
    g_1st_xt_int_dn[x0]=g_lexic2eosub[Index(x0,0,0,0)];
    g_1st_xt_int_up[x0]=g_lexic2eosub[Index(x0,LX-1,0,0)];
    g_1st_xt_ext_dn[x0]=g_lexic2eosub[Index(x0,-1,0,0)];
    g_1st_xt_ext_up[x0]=g_lexic2eosub[Index(x0,LX,0,0)];

    g_1st_yt_int_dn[x0]=g_lexic2eosub[Index(x0,0,0,0)];
    g_1st_yt_int_up[x0]=g_lexic2eosub[Index(x0,0,LY-1,0)];
    g_1st_yt_ext_dn[x0]=g_lexic2eosub[Index(x0,0,-1,0)];
    g_1st_yt_ext_up[x0]=g_lexic2eosub[Index(x0,0,LY,0)];

    g_1st_zt_int_dn[x0]=g_lexic2eosub[Index(x0,0,0,0)];
    g_1st_zt_int_up[x0]=g_lexic2eosub[Index(x0,0,0,LZ-1)];
    g_1st_zt_ext_dn[x0]=g_lexic2eosub[Index(x0,0,0,-1)];
    g_1st_zt_ext_up[x0]=g_lexic2eosub[Index(x0,0,0,LZ)];
  }
# endif
  /* Starting points of the t, x, y, z slices at the borders (eo) */
  g_1st_t_int_dn=g_lexic2eosub[Index(0,0,0,0)];
  g_1st_t_int_up=g_lexic2eosub[Index(T-1,0,0,0)];
  g_1st_t_ext_dn=g_lexic2eosub[Index(-1,0,0,0)];
  g_1st_t_ext_up=g_lexic2eosub[Index(T,0,0,0)];
  
  g_1st_x_int_dn=g_lexic2eosub[Index(0,0,0,0)];
  g_1st_x_int_up=g_lexic2eosub[Index(0,LX-1,0,0)];
  g_1st_x_ext_dn=g_lexic2eosub[Index(0,-1,0,0)];
  g_1st_x_ext_up=g_lexic2eosub[Index(0,LX,0,0)];
  
  g_1st_y_int_dn=g_lexic2eosub[Index(0,0,0,0)];
  g_1st_y_int_up=g_lexic2eosub[Index(0,0,LY-1,0)];
  g_1st_y_ext_dn=g_lexic2eosub[Index(0,0,-1,0)];
  g_1st_y_ext_up=g_lexic2eosub[Index(0,0,LY,0)];
  
  g_1st_z_int_dn=g_lexic2eosub[Index(0,0,0,0)];
  g_1st_z_int_up=g_lexic2eosub[Index(0,0,0,LZ-1)];
  g_1st_z_ext_dn=g_lexic2eosub[Index(0,0,0,-1)];
  g_1st_z_ext_up=g_lexic2eosub[Index(0,0,0,LZ)];

  /* non even-odd */
  gI_0_0_0_0=Index(0,0,0,0);

  gI_L_0_0_0=Index(T,0,0,0);
  gI_Lm1_0_0_0=Index(T-1,0,0,0);
  gI_m1_0_0_0=Index(-1,0,0,0);
  gI_p1_0_0_0=Index(1,0,0,0);
  gI_Lp1_0_0_0=Index(T+1,0,0,0);
  gI_Lm2_0_0_0=Index(T-2,0,0,0);
  gI_m2_0_0_0=Index(-2,0,0,0);

  gI_0_L_0_0=Index(0,LX,0,0);
  gI_0_Lm1_0_0=Index(0,LX-1,0,0);
  gI_0_m1_0_0=Index(0,-1,0,0);
  gI_0_p1_0_0=Index(0,1,0,0);
  gI_0_Lp1_0_0=Index(0,LX+1,0,0);
  gI_0_Lm2_0_0=Index(0,LX-2,0,0);
  gI_0_m2_0_0=Index(0,-2,0,0);

  gI_0_0_L_0=Index(0,0,LY,0);
  gI_0_0_Lm1_0=Index(0,0,LY-1,0);
  gI_0_0_m1_0=Index(0,0,-1,0);
  gI_0_0_p1_0=Index(0,0,1,0);
  gI_0_0_Lp1_0=Index(0,0,LY+1,0);
  gI_0_0_Lm2_0=Index(0,0,LY-2,0);
  gI_0_0_m2_0=Index(0,0,-2,0);

  gI_0_0_0_L=Index(0,0,0,LZ);
  gI_0_0_0_Lm1=Index(0,0,0,LZ-1);
  gI_0_0_0_m1=Index(0,0,0,-1);
  gI_0_0_0_p1=Index(0,0,0,1);
  gI_0_0_0_Lp1=Index(0,0,0,LZ+1);
  gI_0_0_0_Lm2=Index(0,0,0,LZ-2);
  gI_0_0_0_m2=Index(0,0,0,-2);

  gI_L_L_0_0=Index(T,LX,0,0);
  gI_Lm1_L_0_0=Index(T-1,LX,0,0);
  gI_m1_L_0_0=Index(-1,LX,0,0);

  gI_p1_L_0_0=Index(1,LX,0,0);
  gI_Lp1_L_0_0=Index(T+1,LX,0,0);
  gI_Lm2_L_0_0=Index(T-2,LX,0,0);
  gI_m2_L_0_0=Index(-2,LX,0,0);

  gI_L_Lp1_0_0=Index(T,LX+1,0,0);
  gI_Lm1_Lp1_0_0=Index(T-1,LX+1,0,0);
  gI_m1_Lp1_0_0=Index(-1,LX+1,0,0);

  gI_0_L_L_0=Index(0,LX,LY,0);
  gI_0_Lm1_L_0=Index(0,LX-1,LY,0);
  gI_0_m1_L_0=Index(0,-1,LY,0);

  gI_L_0_L_0=Index(T,0,LY,0);
  gI_L_0_Lm1_0=Index(T,0,LY-1,0);
  gI_L_0_m1_0=Index(T,0,-1,0);

  gI_0_p1_L_0=Index(0,1,LY,0);
  gI_0_Lp1_L_0=Index(0,LX+1,LY,0);
  gI_0_Lm2_L_0=Index(0,LX-2,LY,0);
  gI_0_m2_L_0=Index(0,-2,LY,0);

  gI_0_L_Lp1_0=Index(0,LX,LY+1,0);
  gI_0_Lm1_Lp1_0=Index(0,LX-1,LY+1,0);
  gI_0_m1_Lp1_0=Index(0,-1,LY+1,0);

  gI_Lp1_0_L_0=Index(T+1,0,LY,0);
  gI_Lp1_0_Lm1_0=Index(T+1,0,LY-1,0);
  gI_Lp1_0_m1_0=Index(T+1,0,-1,0);

  gI_L_0_p1_0=Index(T,0,1,0);
  gI_L_0_Lp1_0=Index(T,0,LY+1,0);
  gI_L_0_Lm2_0=Index(T,0,LY-2,0);
  gI_L_0_m2_0=Index(T,0,-2,0);

  gI_0_L_0_L=Index(0,LX,0,LZ);
  gI_0_Lm1_0_L=Index(0,LX-1,0,LZ);
  gI_0_m1_0_L=Index(0,-1,0,LZ);

  gI_L_0_0_L=Index(T,0,0,LZ);
  gI_L_0_0_Lm1=Index(T,0,0,LZ-1);
  gI_L_0_0_m1=Index(T,0,0,-1);

  gI_0_L_0_L=Index(0,LX,0,LZ);
  gI_0_Lm1_0_L=Index(0,LX-1,0,LZ);
  gI_0_m1_0_L=Index(0,-1,0,LZ);

  gI_Lp1_0_0_L=Index(T+1,0,0,LZ);
  gI_Lp1_0_0_Lm1=Index(T+1,0,0,LZ-1);
  gI_Lp1_0_0_m1=Index(T+1,0,0,-1);

  gI_L_0_0_p1=Index(T,0,0,1);
  gI_L_0_0_Lp1=Index(T,0,0,LZ+1);
  gI_L_0_0_Lm2=Index(T,0,0,LZ-2);
  gI_L_0_0_m2=Index(T,0,0,-2);

  gI_0_L_0_Lp1=Index(0,LX,0,LZ+1);
  gI_0_Lm1_0_Lp1=Index(0,LX-1,0,LZ+1);
  gI_0_m1_0_Lp1=Index(0,-1,0,LZ+1);

  gI_0_p1_0_L=Index(0,1,0,LZ);
  gI_0_Lp1_0_L=Index(0,LX+1,0,LZ);
  gI_0_Lm2_0_L=Index(0,LX-2,0,LZ);
  gI_0_m2_0_L=Index(0,-2,0,LZ);

  gI_0_0_L_L=Index(0,0,LY,LZ);
  gI_0_0_Lm1_L=Index(0,0,LY-1,LZ);
  gI_0_0_m1_L=Index(0,0,-1,LZ);

  gI_0_0_L_Lp1=Index(0,0,LY,LZ+1);
  gI_0_0_Lm1_Lp1=Index(0,0,LY-1,LZ+1);
  gI_0_0_m1_Lp1=Index(0,0,-1,LZ+1);

  gI_0_0_p1_L=Index(0,0,1,LZ);
  gI_0_0_Lp1_L=Index(0,0,LY+1,LZ);
  gI_0_0_Lm2_L=Index(0,0,LY-2,LZ);
  gI_0_0_m2_L=Index(0,0,-2,LZ);

  gI_Lp1_m1_0_0=Index(T+1,-1,0,0);
  gI_m2_m1_0_0=Index(-2,-1,0,0);
  gI_m2_0_L_0=Index(-2,0,LY,0);
  gI_m2_0_m1_0=Index(-2,0,-1,0);
  gI_0_Lp1_m1_0=Index(0,LX+1,-1,0);
  gI_0_m2_m1_0=Index(0,-2,-1,0);
  gI_m2_0_0_L=Index(-2,0,0,LZ);
  gI_m2_0_0_m1=Index(-2,0,0,-1);
  gI_0_Lp1_0_m1=Index(0,LX+1,0,-1);
  gI_0_m2_0_m1=Index(0,-2,0,-1);
  gI_0_0_Lp1_m1=Index(0,0,LY+1,-1);
  gI_0_0_m2_m1=Index(0,0,-2,-1);
  gI_m1_0_0_m2=Index(-1,0,0,-2);
#endif /* _INDEX_INDEP_GEOM */

#ifdef WITHLAPH
  tempT=T;
  T=1;
  tempV=VOLUME;
  VOLUME=SPACEVOLUME;
  tempR=RAND;
  RAND=SPACERAND;
  gI_0_0_0=Index(0,0,0,0);
  gI_L_0_0=Index(0,LX,0,0);
  gI_Lm1_0_0=Index(0,LX-1,0,0);
  gI_m1_0_0=Index(0,-1,0,0);
  gI_0_L_0=Index(0,0,LY,0);
  gI_0_Lm1_0=Index(0,0,LY-1,0);
  gI_0_m1_0=Index(0,0,-1,0);
  gI_0_0_L=Index(0,0,0,LZ);
  gI_0_0_Lm1=Index(0,0,0,LZ-1);
  gI_0_0_m1=Index(0,0,0,-1);
  T=tempT;
  VOLUME=tempV;
  RAND=tempR;
#endif

#if ( defined PARALLELXYZT || defined PARALLELXYZ )
  check_struct_zt=0;
  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    isp = 0;
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 +  
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) { 
	  g_field_z_ipt_even[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]];
# ifdef _INDEX_INDEP_GEOM
	  g_field_z_disp_even_dn[ix] = g_field_z_ipt_even[ix] - g_1st_z_int_dn;
#  if defined _USE_TSPLITPAR
	  g_field_zt_disp_even_dn[x0][isp] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]]-g_1st_zt_int_dn[x0];
	  if(g_field_zt_disp_even_dn[x0][isp] != g_field_zt_disp_even_dn[x0 % 2][isp]){
	    check_struct_zt=1;
	  }
	  isp++;
#  endif
# endif
	  ix++;
	}
      }
    }
  }
  for(x0 = 0; x0 < T; x0++) {
    isp = 0;
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 + (LZ-1) + 
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==0) { 
	  g_field_z_ipt_even[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]];
# ifdef _INDEX_INDEP_GEOM
	  g_field_z_disp_even_up[ix-T*LX*LY/2] = g_field_z_ipt_even[ix] - g_1st_z_int_up;
#  if defined _USE_TSPLITPAR
	  g_field_zt_disp_even_up[x0][isp] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]]-g_1st_zt_int_up[x0];
	  if(g_field_zt_disp_even_up[x0][isp] != g_field_zt_disp_even_up[x0 % 2][isp]){
	    check_struct_zt=1;
	  }
	  isp++;
#  endif
# endif
	  ix++;
	}
      }
    }
  }
  ix = 0;
  for(x0 = 0; x0 < T; x0++) {
    isp = 0;
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 +  
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==1) { 
	  g_field_z_ipt_odd[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]];
# ifdef _INDEX_INDEP_GEOM
	  g_field_z_disp_odd_dn[ix] = g_field_z_ipt_odd[ix] - g_1st_z_int_dn;
#  if defined _USE_TSPLITPAR
	  g_field_zt_disp_odd_dn[x0][isp] = g_lexic2eosub[ g_ipt[x0][x1][x2][0]]-g_1st_zt_int_dn[x0];
	  if(g_field_zt_disp_odd_dn[x0][isp] != g_field_zt_disp_odd_dn[x0 % 2][isp]){
	    check_struct_zt=1;
	  }
	  isp++;
#  endif
# endif
	  ix++;
	}
      }
    }
  }
  for(x0 = 0; x0 < T; x0++) {
    isp = 0;
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	if((x0 + x1 + x2 + (LZ-1) + 
	    g_proc_coords[0]*T + g_proc_coords[1]*LX +  
	    g_proc_coords[2]*LY + g_proc_coords[3]*LZ)%2==1) { 
	  g_field_z_ipt_odd[ix] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]];
# ifdef _INDEX_INDEP_GEOM
	  g_field_z_disp_odd_up[ix-T*LX*LY/2] = g_field_z_ipt_odd[ix] - g_1st_z_int_up;
#  if defined _USE_TSPLITPAR
	  g_field_zt_disp_odd_up[x0][isp] = g_lexic2eosub[ g_ipt[x0][x1][x2][LZ-1]]-g_1st_zt_int_up[x0];
	  if(g_field_zt_disp_odd_up[x0][isp] != g_field_zt_disp_odd_up[x0 % 2][isp]){
	    check_struct_zt=1;
	  }
	  isp++;
#  endif
# endif
	  ix++;
	}
      }
    }
  }

#  if defined _USE_TSPLITPAR
  if(check_struct_zt !=0){
    if(g_proc_id == 0) {
      fprintf(stderr,"Error in assuming the structure of dispacements of zt slice\n");
      fflush(stderr);
      MPI_Finalize();
      exit(-1);
    }
  }
#  endif

# ifdef _INDEX_INDEP_GEOM
  /* Define the MPI_Types using the displacement vectors above and MPI_Type_indexed() */
  ones=malloc(T*LX*LY/2*sizeof(int));
  for(j=0;j<T*LX*LY/2;j++) ones[j]=1;
  MPI_Type_indexed(T*LX*LY/2,ones,g_field_z_disp_even_dn,field_point,&field_z_slice_even_dn);
  MPI_Type_indexed(T*LX*LY/2,ones,g_field_z_disp_even_up,field_point,&field_z_slice_even_up);
  MPI_Type_indexed(T*LX*LY/2,ones,g_field_z_disp_odd_dn,field_point,&field_z_slice_odd_dn);
  MPI_Type_indexed(T*LX*LY/2,ones,g_field_z_disp_odd_up,field_point,&field_z_slice_odd_up);
  MPI_Type_commit(&field_z_slice_even_dn);
  MPI_Type_commit(&field_z_slice_even_up);
  MPI_Type_commit(&field_z_slice_odd_dn);
  MPI_Type_commit(&field_z_slice_odd_up);
  free(ones);

#  if defined _USE_TSPLITPAR
  /* LZ and T*LX*LY are required to be even, but LX*LY does not need to */
  /* length zt slice=ceil(LX*LY/2) = (LX*LY+1)/2, if parity(t)*parity(z)*globalparity=1   */
  /* length zt slice=floor(LX*LY/2) = LX*LY/2, if parity(t)*parity(z)*globalparity=-1   */
  lsliceL=(LX*LY+1)/2;
  lsliceS=(LX*LY/2);
  oneL=malloc(lsliceL*sizeof(int));
  oneS=malloc(lsliceS*sizeof(int));
  for(j=0;j<lsliceL;j++) oneL[j]=1;
  for(j=0;j<lsliceS;j++) oneS[j]=1;

  MPI_Type_indexed(lsliceL,oneL,g_field_zt_disp_even_dn[0],field_point,&field_zt_slice_even_dn_et);
  MPI_Type_commit(&field_zt_slice_even_dn_et);
  MPI_Type_indexed(lsliceS,oneS,g_field_zt_disp_even_up[0],field_point,&field_zt_slice_even_up_et);
  MPI_Type_commit(&field_zt_slice_even_up_et);
  MPI_Type_indexed(lsliceS,oneS,g_field_zt_disp_odd_dn[0],field_point,&field_zt_slice_odd_dn_et);
  MPI_Type_commit(&field_zt_slice_odd_dn_et);
  MPI_Type_indexed(lsliceL,oneL,g_field_zt_disp_odd_up[0],field_point,&field_zt_slice_odd_up_et);
  MPI_Type_commit(&field_zt_slice_odd_up_et);
  MPI_Type_indexed(lsliceS,oneS,g_field_zt_disp_even_dn[1],field_point,&field_zt_slice_even_dn_ot);
  MPI_Type_commit(&field_zt_slice_even_dn_ot);
  MPI_Type_indexed(lsliceL,oneL,g_field_zt_disp_even_up[1],field_point,&field_zt_slice_even_up_ot);
  MPI_Type_commit(&field_zt_slice_even_up_ot);
  MPI_Type_indexed(lsliceL,oneL,g_field_zt_disp_odd_dn[1],field_point,&field_zt_slice_odd_dn_ot);
  MPI_Type_commit(&field_zt_slice_odd_dn_ot);
  MPI_Type_indexed(lsliceS,oneS,g_field_zt_disp_odd_up[1],field_point,&field_zt_slice_odd_up_ot);
  MPI_Type_commit(&field_zt_slice_odd_up_ot);

#  endif
# endif

#endif /* PARALLELXYZ || PARALLELXYZT*/

  /* The rectangular gauge action part */
  /* Everything is stored behind VOLUMEPLUSRAND-1 !*/
  if(g_dbw2rand != 0) {
    if(g_proc_id == 0) {
      printf("# Initialising rectangular gauge action stuff\n");
      fflush(stdout);
    }
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* t2 Rand and t2x and t2y */
	    x0 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### -2t %d %d %d %d\n",x0, x1, x2, x3);
	    }

	    g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    g_idn[ix][0] = -1;

	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	    x0 = T+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### +2t %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    g_iup[ix][0] = -1;
	    g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	  }
	}
      }
    }    
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* x2-Rand and x2t and x2y */
	    x1 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### -2x %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    g_idn[ix][1] = -1;

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);

	    x1 = LX+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND) {
	      printf("#### +2x %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);

	    g_iup[ix][1] = -1;
	    g_idn[ix][1] = Index(x0, x1-1, x2, x3);

	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);

	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++) {
	for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x3 < 0 || x3 > LZ-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* y2-Rand y2t and y2x */
	    x2 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### -2y %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    g_idn[ix][2] = -1;
	    
	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	    
	    x2 = LY+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### +2y %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 < T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    g_iup[ix][2] = -1;
	    g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    if(x3 < LZ) g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    if(x3 > -1) g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
#if (defined PARALLELXYZT || defined PARALLELXYZ)
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++) {
	for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++) {
	  bndcnt = 0;
	  if(x0 < 0 || x0 > T-1) bndcnt++;
	  if(x1 < 0 || x1 > LX-1) bndcnt++;
	  if(x2 < 0 || x2 > LY-1) bndcnt++;	  
	  if(bndcnt < 2) {
	    /* z2-Rand t2z and z2x z2y*/
	    x3 = -2;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### -2z %d %d %d %d %d %d %d %d %d %d %d\n",x0, x1, x2, x3, ix, 
		     VOLUMEPLUSRAND, VOLUMEPLUSRAND + g_dbw2rand, T, LX, LY, LZ);
	    }
	    if(x0 <  T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    g_iup[ix][3] = Index(x0, x1, x2, x3+1);
	    g_idn[ix][3] = -1;
	    
	    x3 = LZ+1;
	    ix = Index(x0, x1, x2, x3);
	    if(ix < VOLUMEPLUSRAND || ix >= VOLUMEPLUSRAND + g_dbw2rand) {
	      printf("#### +2z %d %d %d %d\n",x0, x1, x2, x3);
	    }
	    if(x0 <  T) g_iup[ix][0] = Index(x0+1, x1, x2, x3);
	    if(x0 > -1) g_idn[ix][0] = Index(x0-1, x1, x2, x3);
	    
	    if(x1 < LX) g_iup[ix][1] = Index(x0, x1+1, x2, x3);
	    if(x1 > -1) g_idn[ix][1] = Index(x0, x1-1, x2, x3);
	    
	    if(x2 < LY) g_iup[ix][2] = Index(x0, x1, x2+1, x3);
	    if(x2 > -1) g_idn[ix][2] = Index(x0, x1, x2-1, x3);
	    
	    g_iup[ix][3] = -1;
	    g_idn[ix][3] = Index(x0, x1, x2, x3-1);
	  }
	}
      }
    }
#endif
  }

  Hopping_Matrix_Indices();

  free(xeven);
}


void Hopping_Matrix_Indices(){
  int ix;
  int ioff = (VOLUME+RAND)/2;
  /**************** loop over all lattice sites ****************/
  for (int icx = 0, icy = (VOLUME+RAND)/2; icx < VOLUME/2; icx++, icy++)
  {
    ix=g_eo2lexic[icx];
    /*********************** direction +0 ************************/
    g_hi[(16*icx)] = g_iup[ix][0];
    g_hi[(16*icx)+1] = g_lexic2eosub[g_hi[(16*icx)]];
    g_hi[(16*icx)] = ix;
    /*********************** direction -0 ************************/
    g_hi[(16*icx)+2] = g_idn[ix][0];
    g_hi[(16*icx)+3] = g_lexic2eosub[g_hi[(16*icx)+2]];
    /*********************** direction +1 ************************/
    g_hi[(16*icx)+4] = g_iup[ix][1];
    g_hi[(16*icx)+5] = g_lexic2eosub[g_hi[(16*icx)+4]];
    /*********************** direction -1 ************************/
    g_hi[(16*icx)+6] = g_idn[ix][1];
    g_hi[(16*icx)+7] = g_lexic2eosub[g_hi[(16*icx)+6]];
    /*********************** direction +2 ************************/
    g_hi[(16*icx)+8] = g_iup[ix][2];
    g_hi[(16*icx)+9] = g_lexic2eosub[g_hi[(16*icx)+8]];
    /*********************** direction -2 ************************/
    g_hi[(16*icx)+10] = g_idn[ix][2];
    g_hi[(16*icx)+11] = g_lexic2eosub[g_hi[(16*icx)+10]];
    /*********************** direction +3 ************************/
    g_hi[(16*icx)+12] = g_iup[ix][3];
    g_hi[(16*icx)+13] = g_lexic2eosub[g_hi[(16*icx)+12]];
    /*********************** direction -3 ************************/
    g_hi[(16*icx)+14] = g_idn[ix][3];
    g_hi[(16*icx)+15] = g_lexic2eosub[g_hi[(16*icx)+14]];
    /************************ end of loop ************************/
    ix=g_eo2lexic[icx+ioff];
    /*********************** direction +0 ************************/
    g_hi[(16*icy)] = g_iup[ix][0];
    g_hi[(16*icy)+1] = g_lexic2eosub[g_hi[(16*icy)]];
    g_hi[(16*icy)] = ix;
    /*********************** direction -0 ************************/
    g_hi[(16*icy)+2] = g_idn[ix][0];
    g_hi[(16*icy)+3] = g_lexic2eosub[g_hi[(16*icy)+2]];
    /*********************** direction +1 ************************/
    g_hi[(16*icy)+4] = g_iup[ix][1];
    g_hi[(16*icy)+5] = g_lexic2eosub[g_hi[(16*icy)+4]];
    /*********************** direction -1 ************************/
    g_hi[(16*icy)+6] = g_idn[ix][1];
    g_hi[(16*icy)+7] = g_lexic2eosub[g_hi[(16*icy)+6]];
    /*********************** direction +2 ************************/
    g_hi[(16*icy)+8] = g_iup[ix][2];
    g_hi[(16*icy)+9] = g_lexic2eosub[g_hi[(16*icy)+8]];
    /*********************** direction -2 ************************/
    g_hi[(16*icy)+10] = g_idn[ix][2];
    g_hi[(16*icy)+11] = g_lexic2eosub[g_hi[(16*icy)+10]];
    /*********************** direction +3 ************************/
    g_hi[(16*icy)+12] = g_iup[ix][3];
    g_hi[(16*icy)+13] = g_lexic2eosub[g_hi[(16*icy)+12]];
    /*********************** direction -3 ************************/
    g_hi[(16*icy)+14] = g_idn[ix][3];
    g_hi[(16*icy)+15] = g_lexic2eosub[g_hi[(16*icy)+14]];
    /************************ end of loop ************************/

  }
  g_hi[(16*(VOLUME+RAND))] = 0;
  g_hi[(16*(VOLUME+RAND))+1] = 0;
  return;
}

