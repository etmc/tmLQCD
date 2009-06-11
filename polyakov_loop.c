/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2005 Urs Wenger
 *               2008 Marcus Petschlies
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
 * Routine to calculate the Polyakov loop.
 * 
 * Author: Urs Wenger <urs.wenger@desy.de>
 * Date: January 2005
 *
 * Polyakov loop in time direction added by Marcus Petschlies
 *       2008
 *
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "sse.h"
#include "su3.h"
#include "mpi_init.h"
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
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 10);
    MPI_Finalize();
#endif
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


/* here comes the one in time direction */

int polyakov_loop_0(const int nstore, complex *pl) {

  int i0, i1, i2, i3, ixyz, ixyzt, ixyzt_up, VOL3, VOLUME3;
  int L0, L1, L2, L3;
  complex pl_tmp, tr, ts, tt, kc, ks;
  su3 *tmp_loc = NULL, tmp, tmp2;
  su3 *v = NULL, *w = NULL;
  
  FILE *ofs = NULL;
  
#ifdef MPI
  int iproc;
  MPI_Status status;
  su3 *tmp_nnb = NULL;
#endif
  
  L0 = LX; /* enable transparent comparison with existing Polyakov routines */
  L1 = LY; /* in spatial directions                                         */
  L2 = LZ;
  L3 = T;
  
  /**************
   * local part *
   **************/
  VOL3 = L0*L1*L2;
  tmp_loc = (su3 *)calloc(VOL3, sizeof(su3));
  
  for(i0 = 0; i0 < LX; i0++) {
    for(i1 = 0; i1 < LY; i1++) {
      for(i2 = 0; i2 < LZ; i2++) {
	ixyz = (i2 * L1 + i1) * L0 + i0;
	i3 = 0;
	ixyzt    = g_ipt[i3][i0][i1][i2];
	ixyzt_up = g_iup[ixyzt][0];
	v = &g_gauge_field[ixyzt][0];
	w = &g_gauge_field[ixyzt_up][0];
	_su3_times_su3(tmp, *v, *w);
	
	for(i3 = 1; i3 < L3-1; i3++) {
	  _su3_assign(tmp2,tmp);
	  ixyzt_up = g_iup[ixyzt_up][0];
	  w = &g_gauge_field[ixyzt_up][0];
	  _su3_times_su3(tmp, tmp2, *w);
	}
	_su3_assign(tmp_loc[ixyz],tmp);
      }
    }
  }
  
  /********************************************************************************/
    
#ifdef MPI
  /***************
   * global part *
   ***************/
  /* (1) collect contributions from different time slices to nodes with t-coord. 0 */
  tmp_nnb = (su3*)calloc(VOL3, sizeof(su3)); /* contains the next-neighbour-part*/

  /* note: in the following loop t is taken as the time coordinate of nodes */
  for(iproc = g_nproc_t-1; iproc > 0; iproc--) {
    if(g_proc_coords[0] == iproc) /* node is in the {t=iproc}-hyperplane */ {
      MPI_Send(tmp_loc, VOL3, mpi_su3, g_nb_t_dn, 100+g_cart_id, g_cart_grid);
      /* send tmp_loc from {t=iproc}-hyperplane to {t=iproc-1}-hyperplane */
    }
    if(g_proc_coords[0] == iproc-1) {
      /* so the node is right below the sending one in time(= 0)-direction */
      MPI_Recv(tmp_nnb, VOL3, mpi_su3, g_nb_t_up, 100+g_nb_t_up, g_cart_grid, &status); 
      /* receive tmp_loc from the tmp_loc from the
	 {t=my_own_t_index+1}-hyperplane */
      for(ixyz=0; ixyz<VOL3; ixyz++) {
	/* multiply all matrices in tmp_nbb to my own in tmp_loc from the right */
	v = tmp_loc+ixyz; 
	w = tmp_nnb+ixyz; 
	_su3_assign(tmp2, *v);
	_su3_times_su3(*v, tmp2, *w);
      }
    }
    /* if iproc==0 then the node with g_proc_coords[0]=0 will finally contain
       the product of all contributions from all {t=const.}-planes */
  }
  
  /* (2) nodes with time coordinate 0 sum traces over local spatial points */
#endif
  _complex_zero(pl_tmp);
  /*   pl_tmp.re = 0.0; pl_tmp.im = 0.0; */
  if(g_proc_coords[0] == 0) {
    
    kc.re = 0.0; kc.im = 0.0; ks.re = 0.0; ks.im = 0.0;
    for(ixyz = 0; ixyz < VOL3; ixyz++) /* Kahan-summation of traces */ {
      pl_tmp.re = (tmp_loc[ixyz]).c00.re + (tmp_loc[ixyz]).c11.re + (tmp_loc[ixyz]).c22.re;
      pl_tmp.im = (tmp_loc[ixyz]).c00.im + (tmp_loc[ixyz]).c11.im + (tmp_loc[ixyz]).c22.im;
      tr.re=pl_tmp.re+kc.re;
      ts.re=tr.re+ks.re;
      tt.re=ts.re-ks.re;
      ks.re=ts.re;
      kc.re=tr.re-tt.re;
      tr.im=pl_tmp.im+kc.im;
      ts.im=tr.im+ks.im;
      tt.im=ts.im-ks.im;
      ks.im=ts.im;
      kc.im=tr.im-tt.im;
    }
    pl_tmp.re = ks.re + kc.re;
    pl_tmp.im = ks.im + kc.im;
  }
  
#ifdef MPI
  /* (3) sum over all contributions from all nodes (also nodes with pl_tmp=0;
     apparently the easiest way) */
  MPI_Reduce(&pl_tmp, pl, 1, MPI_COMPLEX, MPI_SUM, 0, g_cart_grid);
  /*   MPI_Reduce(&(pl_tmp.re), &((*pl).re), 1, MPI_DOUBLE, MPI_SUM, 0, g_cart_grid);  */
  /*   MPI_Reduce(&(pl_tmp.im), &((*pl).im), 1, MPI_DOUBLE, MPI_SUM, 0, g_cart_grid);  */
#else
  (*pl).re = pl_tmp.re;
  (*pl).im = pl_tmp.im;
#endif
  
  /* normalization */
  VOLUME3 = VOL3;

  if(g_proc_id == 0) {
    VOLUME3 = VOLUME3 * g_nproc_x*g_nproc_y*g_nproc_z;
    (*pl).re /= 3*VOLUME3;
    (*pl).im /= 3*VOLUME3;
  }
  
  /* write result to file */
  if (g_proc_id == 0) {
    if (nstore == 0) {
      ofs = fopen("polyakov_loop_0.dat","w");
    }
    else {
      ofs = fopen("polyakov_loop_0.dat","a");
    }
    fprintf(ofs, "%12.10e %12.10e\n", (*pl).re, (*pl).im); 
    fclose(ofs);
  }
#ifdef MPI
  free(tmp_nnb);
#endif
  free(tmp_loc);
  return(0);
}
