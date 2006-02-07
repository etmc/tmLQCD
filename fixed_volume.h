/* $Id$ */
/********************************************************
 *
 * In case the code was configure for fixed volume at
 * compiletime, the parameters have to be set here!
 *
 *
 ********************************************************/

#ifndef _FIXED_VOLUME_H
#define _FIXED_VOLUME_H

#  if defined FIXEDVOLUME

/* Set the next 8 number! */

#    define Tdef 48
#    define Xdef 24
#    define Ydef 24
#    define Zdef 24

#    define N_PROC_T 1
#    define N_PROC_X 1
#    define N_PROC_Y 1
#    define N_PROC_Z 1

/* The rest is done automatially */

#    define T (Tdef/N_PROC_T)
#    define LX (Xdef/N_PROC_X)
#    define LY (Ydef/N_PROC_Y)
#    define LZ (Zdef/N_PROC_Z)
#    define VOLUME (T*LX*LY*LZ)

#    ifdef PARALLELT  
#      define RAND (2*LX*LY*LZ)
#      define EDGES 0
#    elif defined PARALLELXT
#      define RAND (2*LZ*(LY*LX + T*LY))
#      define EDGES (4*LZ*LY);
#    elif defined PARALLELXYT
#      define RAND (2*LZ*(LY*LX + T*LY + T*LX))
#      define EDGES (4*LZ*(LY + T + LX))
#    elif defined PARALLELXYZT
#      define RAND (2*LZ*(LY*LX + T*LY + T*LX) + 2*T*LX*LY)
#      define EDGES (4*LZ*(LY + T + LX) + 4*LY*T + 4*LY*LX + 4*T*LX)
#    else
#      define RAND 0
#      define EDGES 0
#    endif
  /* Note that VOLUMEPLUSRAND is in general not equal to VOLUME+RAND */
  /* VOLUMEPLUSRAND rather includes the edges */
#    define VOLUMEPLUSRAND (VOLUME + RAND + EDGES)

#  endif

#endif
