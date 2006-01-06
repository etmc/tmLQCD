/* $Id$ */

/*******************************************************************
 *
 * Routine to calculate the Polyakov loop.
 * 
 * Author: Urs Wenger <urs.wenger@desy.de>
 * Date: January 2005
 *
*******************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "polyakov_loop.h"

void polyakov_loop(complex * pl_, const int mu) {

  static int i0, i1, i2, i3, L0, L1, L2, L3, ixyzt, ixyzt_up;
  static double vol;
  static su3 tmp, tmp2; 
  su3 *v = NULL , *w = NULL;
  static complex pl; 
  /* For the Kahan summation:*/
#ifdef MPI
  static complex pls; 
#endif
  static complex ks, kc, tr, ts, tt;
  kc.re=0.0; ks.re=0.0;
  kc.im=0.0; ks.im=0.0;
  

  /* For the moment only the Polyakov loop in y- and z-direction 
     are implemented, since they are not affected by parallelisation: */
  if(mu == 0 || mu == 1 || mu > 3) {
      fprintf(stderr, "Wrong parameter for Polyakov loop calculation in polyakov_loop.c:\n");
      fprintf(stderr, "Only direction %d and %d are allowed.\n",2,3);
      fprintf(stderr, "Actual value is %d! Aborting...\n",mu);
      exit(0);
  }


  L0=T;
  L1=LX;
  if(mu==2) {
      L2=LZ;
      L3=LY;
  }
  else {
      L2=LY;
      L3=LZ;
  }
  /* loop over the spatial sites: */
  for (i0=0; i0 < L0; i0++) {
    for (i1=0; i1 < L1; i1++) {
      for (i2=0; i2 < L2; i2++) {
	/* at each spatial site multiply the links in 
	   temporal direction: */
	  i3 = 0;
	  /* get the site index: */
	  if(mu==2) {
	      ixyzt = g_ipt[i0][i1][i3][i2];
	  }
	  else {
	      ixyzt = g_ipt[i0][i1][i2][i3];
	  }
	  /* and its neigbour in direction mu: */
	  ixyzt_up = g_iup[ixyzt][mu];
	  
	  /* Get the links and multiply them: ixyzt --> ixyzt_up --> */
	  v = &g_gauge_field[ixyzt][mu];
	  w = &g_gauge_field[ixyzt_up][mu];
	  _su3_times_su3(tmp, *v, *w);

	  /* now start the loop over indices in mu-direction: */
	  for (i3=1; i3 < L3-2; i3++) {
	    /* store the current result in v:*/
	    _su3_assign(tmp2,tmp);
	    /* get the next site index: */
	    ixyzt_up = g_iup[ixyzt_up][mu];
	    /* and the corresponding link matrix: */
	    w = &g_gauge_field[ixyzt_up][mu];
	    /* and multiply them: */
	    _su3_times_su3(tmp, tmp2, *w);
	  }
	  
	  /* for the last link we directly take the complex trace: */
	  ixyzt_up = g_iup[ixyzt_up][mu];
	  w = &g_gauge_field[ixyzt_up][mu];
	  _trace_su3_times_su3(pl,tmp,*w);

	  /* printf("i0=%d, i1=%d, i2=%d, pl=(%e,%e)\n",i0,i1,i2,pl.re,pl.im);*/

	  /* Kahan summation for real and imaginary part: */
	  tr.re=pl.re+kc.re;
	  ts.re=tr.re+ks.re;
	  tt.re=ts.re-ks.re;
	  ks.re=ts.re;
	  kc.im=tr.im-tt.im;
	  tr.im=pl.im+kc.im;
	  ts.im=tr.im+ks.im;
	  tt.im=ts.im-ks.im;
	  ks.im=ts.im;
	  kc.im=tr.im-tt.im;
      }
    }
  }
  /* Finish Kahan summation: */
  /* (Division by 3 is for normalising the colour trace.) */
  pl.re=(kc.re+ks.re)/3.0;
  pl.im=(kc.im+ks.im)/3.0;
/*  printf("Polyakov loop before normalisation, pl.re=%e, pl.im=%e\n",pl.re,pl.im);*/


  /* Collect the results and return:*/
#ifdef MPI
  MPI_Allreduce(&pl, &pls, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  pl=pls;
#endif

  /* Normalise, i.e. divide by the number of loops: */
  vol = (double) L0*L1*L2*g_nproc_t*g_nproc_x;
/*    printf("L0*L1*L2=%d, vol=%e\n",L0*L1*L2,vol);  */
  _div_real(pl,pl,vol);
/*    printf("Polyakov loop after normalisation, pl.re=%e, pl.im=%e\n",pl.re,pl.im) */; 
/*   return pl; */
  (*pl_).re = pl.re;
  (*pl_).im = pl.im;
}
