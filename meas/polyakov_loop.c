/***********************************************************************
 *
 * Copyright (C) 2005 Urs Wenger
 *               2008,2009 Marcus Petschlies
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
#include <time.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <complex.h>
#include "sse.h"
#include "su3.h"
#include "read_input.h"
#include "start.h"
#include "mpi_init.h"
#include "polyakov_loop.h"
#include "gettime.h"

void polyakov_loop(_Complex double * pl_, const int mu) {

  static int i0, i1, i2, i3, L0, L1, L2, L3, ixyzt, ixyzt_up;
  static double vol;
  static su3 tmp, tmp2; 
  su3 *v = NULL , *w = NULL;
  static _Complex double pl; 
  /* For the Kahan summation:*/
#ifdef MPI
  static _Complex double pls; 
#endif
  static _Complex double ks = 0.0, kc = 0.0, tr, ts, tt;
  
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
	
	/* for the last link we directly take the _Complex double trace: */
	ixyzt_up = g_iup[ixyzt_up][mu];
	w = &g_gauge_field[ixyzt_up][mu];
	_trace_su3_times_su3(pl,tmp,*w);
	
	/* printf("i0=%d, i1=%d, i2=%d, pl=(%e,%e)\n",i0,i1,i2,creal(pl),cimag(pl));*/
	
	/* Kahan summation for real and imaginary part: */
	tr = pl + kc;
	ts = tr + ks;
	tt = ts - ks;
	ks = ts;
	kc = tr - tt;
      }
    }
  }
  /* Finish Kahan summation: */
  /* (Division by 3 is for normalising the colour trace.) */
  pl = (kc + ks) / 3.0;
  /*  printf("Polyakov loop before normalisation, pl.re=%e, pl.im=%e\n",creal(pl),cimag(pl));*/
  
  
  /* Collect the results and return:*/
#ifdef MPI
  MPI_Allreduce(&pl, &pls, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  pl=pls;
#endif
  
  /* Normalise, i.e. divide by the number of loops: */
  vol = (double) L0*L1*L2*g_nproc_t*g_nproc_x;
  /*    printf("L0*L1*L2=%d, vol=%e\n",L0*L1*L2,vol);  */
  pl /= vol;
  /*    printf("Polyakov loop after normalisation, pl.re=%e, pl.im=%e\n",creal(pl),cimag(pl)) */; 
  /*   return pl; */
  *pl_ = pl;
}


/* here comes the one in time direction */

int polyakov_loop_0(const int nstore, _Complex double *pl) {

  int i0, i1, i2, i3, ixyz, ixyzt, ixyzt_up, VOL3, VOLUME3;
  int L0, L1, L2, L3;
  double retime, ratime;
  _Complex double pl_tmp, tr, ts, tt, kc, ks;
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
    ratime = gettime();

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
  retime = gettime();
  if(g_debug_level>0) {
    fprintf(stdout, "[polyakov_loop_0 | %3d] time for calculating local part = %e seconds\n", g_cart_id, retime-ratime);
  }
  
  /********************************************************************************/
    
#ifdef MPI
  /***************
   * global part *
   ***************/

  ratime = MPI_Wtime();

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
 
  retime = MPI_Wtime();
  if(g_proc_id==0 && g_debug_level>0) {
    fprintf(stdout, "[polyakov_loop_0 | %3d] time for calculating global part = %e seconds\n", g_cart_id, retime-ratime);
  }
 
  /* (2) nodes with time coordinate 0 sum traces over local spatial points */
#endif
  pl_tmp = 0.0;
  if(g_proc_coords[0] == 0) {
    
    kc = 0.0; ks = 0.0;
    for(ixyz = 0; ixyz < VOL3; ixyz++) /* Kahan-summation of traces */
    {
      pl_tmp = tmp_loc[ixyz].c00 + tmp_loc[ixyz].c11 + tmp_loc[ixyz].c22;
      tr = pl_tmp + kc;
      ts = tr + ks;
      tt = ts - ks;
      ks = ts;
      kc = tr - tt;
    }
    pl_tmp = ks + kc;
  }
  
#ifdef MPI
  /* (3) sum over all contributions from all nodes (also nodes with pl_tmp=0;
     apparently the easiest way) */
  MPI_Reduce(&pl_tmp, pl, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, g_cart_grid);
  /*   MPI_Reduce(&(creal(pl_tmp)), &(pl->re), 1, MPI_DOUBLE, MPI_SUM, 0, g_cart_grid);  */
  /*   MPI_Reduce(&(cimag(pl_tmp)), &(pl->im), 1, MPI_DOUBLE, MPI_SUM, 0, g_cart_grid);  */
#else
  *pl = pl_tmp;
#endif
  
  /* normalization */
  VOLUME3 = VOL3;

  if(g_proc_id == 0)
  {
    VOLUME3 = VOLUME3 * g_nproc_x*g_nproc_y*g_nproc_z;
    *pl /= 3 * VOLUME3;
  }
  
  /* write result to file */
  if (g_proc_id == 0) {
    if (nstore == 0) {
      ofs = fopen("polyakov_loop_0.dat","w");
    }
    else {
      ofs = fopen("polyakov_loop_0.dat","a");
    }
    fprintf(ofs, "%25.16e\t%25.16e\n", creal(*pl), cimag(*pl)); 
    fclose(ofs);
  }
#ifdef MPI
  free(tmp_nnb);
#endif
  free(tmp_loc);
  return(0);
}


/*********************************************************************************/

/* here comes the version using reduction operations for time- (dir==0) or
   z- (dir==3) direction
   the reduction operation is defined in mpi_init.h
*/
void polyakov_loop_measurement(const int nstore, const int id, const int ieo) {
  polyakov_loop_dir(nstore, measurement_list[id].direction);
}


int polyakov_loop_dir(
		      const int nstore /* in  */,
		      const int dir    /* in  */) {

  int ixyz, ixyzt, ixyzt_up, VOL3, VOLUME3, ix, iy, iz, it;
  _Complex double pl_tmp, tr, ts, tt, kc, ks, pl;
  su3 *tmp_loc, tmp, tmp2;
  su3 *u, *v, *w;
  double ratime, retime;
  char filename[50];

  FILE *ofs;

#ifdef MPI
  int rank_slice, rank_ray;
  MPI_Comm slice, ray;
  su3 *tmp_ray;
#endif

  if(dir!=0 && dir!=3 && g_proc_id==0) {
    fprintf(stderr, "Wrong direction; must be 0 (t) or 3 (z)\n");
    return(-1);
  }

  pl = 0.0;

  /********************************************************************************/

  /**************
   * local part *
   **************/
  ratime = gettime();

  if(dir==0) {
    VOL3 = LX*LY*LZ;
    tmp_loc = (su3 *)calloc(VOL3, sizeof(su3));
    if((void*)tmp_loc == NULL) {
      fprintf(stderr, "[%2d] Could not allocate memory for tmp_loc\n", g_proc_id);
      return(-1);
    }

    for(ix=0; ix<LX; ix++) {
      for(iy=0; iy<LY; iy++) {
	for(iz=0; iz<LZ; iz++) {
	  /* ixyz = ix*LY*LZ + iy*LZ + iz */
	  ixyz = (ix * LY + iy) * LZ + iz;
	  it = 0;
	  ixyzt    = g_ipt[it][ix][iy][iz];
	  ixyzt_up = g_iup[ixyzt][0];
	  v = &g_gauge_field[ixyzt][0];
	  w = &g_gauge_field[ixyzt_up][0];
	  u = &tmp;
	  _su3_times_su3(*u, *v, *w);
	  v = &tmp2;
	  for(it=1; it<T-2; it++) {
	    /* swap u and v via w */
	    w = u; u = v; v = w;
	    ixyzt_up = g_iup[ixyzt_up][0];
	    w = &g_gauge_field[ixyzt_up][0];
	    _su3_times_su3(*u, *v, *w);
	  }
	  /* last multiplication for it=T-1 */
	  ixyzt_up = g_iup[ixyzt_up][0];
	  w = &g_gauge_field[ixyzt_up][0];
	  _su3_times_su3(tmp_loc[ixyz],*u, *w);
	}
      }
    }
  }
  else { /* z-direction <=> dir==3 */
    VOL3 = T*LX*LY;
    tmp_loc = (su3 *)calloc(VOL3, sizeof(su3));
    if((void*)tmp_loc == NULL) {
      /* Abort */
    }

    for(it=0; it<T;  it++) {
      for(ix=0; ix<LX; ix++) {
	for(iy=0; iy<LY; iy++) {
	  /* ixyz = it*LX*LY + ix*LY + iy */
	  ixyz = (it * LX + ix) * LY + iy;
	  iz = 0;
	  ixyzt    = g_ipt[it][ix][iy][iz];
	  ixyzt_up = g_iup[ixyzt][3];
	  v = &g_gauge_field[ixyzt][3];
	  w = &g_gauge_field[ixyzt_up][3];
	  u = &tmp;
	  _su3_times_su3(*u, *v, *w);
	  v = &tmp2;
	  for(iz=1; iz<LZ-2; iz++) {
	    /* swap u and v via w */
	    w = u; u = v; v = w;
	    ixyzt_up = g_iup[ixyzt_up][3];
	    w = &g_gauge_field[ixyzt_up][3];
	    _su3_times_su3(*u, *v, *w);
	  }
	  ixyzt_up = g_iup[ixyzt_up][3];
	  w = &g_gauge_field[ixyzt_up][3];
	  _su3_times_su3(tmp_loc[ixyz], *u, *w);
	}
      }
    }

  }
  retime = gettime();
  if(g_debug_level > 0 && g_proc_id == 0) {
    fprintf(stdout, "# [pl02 dir%1d proc%.2d] time for calculating local part"\
	    " = %e seconds\n", dir, g_cart_id, retime-ratime);
  }

  /********************************************************************************/

#ifdef MPI
  /***************
   * global part *
   ***************/
  /* choose the slice and ray communicators according to direction */
  if(dir==0) {
    slice      = g_mpi_time_slices;
    ray        = g_mpi_SV_slices;
    rank_slice = g_mpi_time_rank;
    rank_ray   = g_mpi_SV_rank;
  }
  else {
    slice      = g_mpi_z_slices;
    ray        = g_mpi_ST_slices;
    rank_slice = g_mpi_z_rank;
    rank_ray   = g_mpi_ST_rank;
  }
 
  ratime = MPI_Wtime();

  /* (1) collect contributions from different time/z slices to nodes with rank=0 
     in spatial volume/space-time slices */
#  ifndef PARALLELXYZT
  if(dir==0) {
#  endif
    tmp_ray = (su3*)calloc(VOL3, sizeof(su3)); /* */
    if((void*)tmp_ray== NULL) {
      fprintf(stderr, "[%2d] Could not allocate memory for tmp_ray\n", g_proc_id);
      return(-1);
    }

    MPI_Reduce(tmp_loc, tmp_ray, VOL3, mpi_su3, mpi_reduce_su3_ray, 0, ray);
#  ifndef PARALLELXYZT
  }
#  endif


  retime = MPI_Wtime();
  if(g_proc_id==0 && g_debug_level>0) {
    fprintf(stdout, "# [pl02 dir%1d proc%.2d] time for calculating global part"\
	    " = %e seconds\n", dir, g_cart_id, retime-ratime);
  }

  if(rank_ray == 0) {

#endif
    pl_tmp = 0.0;
    kc = 0.0;
    ks = 0.0;

#ifdef MPI
#  ifdef PARALLELXYZT
    u = tmp_ray;
#  else
    if(dir==0) { u = tmp_ray; }
    else       { u = tmp_loc; }
#  endif
#else
    u = tmp_loc;
#endif

    for(ixyz=0; ixyz<VOL3; ixyz++) /* Kahan-summation of traces */
    {
      pl_tmp = u[ixyz].c00 + u[ixyz].c11 + u[ixyz].c22;
      tr = pl_tmp + kc;
      ts = tr + ks;
      tt = ts - ks;
      ks = ts;
      kc = tr - tt;
    }
    pl_tmp = ks + kc;

#ifdef MPI
    MPI_Reduce(&pl_tmp, &pl, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, slice);
  }
#  ifndef PARALLELXYZT
  if(dir==0) {
#  endif
    free(tmp_ray);
#  ifndef PARALLELXYZT
  }
#  endif
    
#else
  pl = pl_tmp;
#endif

  /* normalization pl |-> pl / ( 3 * 3-dim. volume)*/
  VOLUME3 = VOL3;
   
#ifdef MPI
  if(rank_slice==0 && rank_ray==0) { /* this process has the sum 
					of the Polyakov loop values */
    if(dir==0) { 
      VOLUME3 = VOLUME3 * g_nproc_x*g_nproc_y*g_nproc_z;
    }
    else {
      VOLUME3 = VOLUME3 * g_nproc_t*g_nproc_x*g_nproc_y;
    }
#endif
    pl /= 3. * VOLUME3;

    /* write result to file */
    sprintf(filename, "polyakovloop_dir%1d", dir);
    if (nstore == 0) {
      ofs = fopen(filename,"w");
    }
    else {
      ofs = fopen(filename,"a");
    }
    if((void*)ofs == NULL) {
      fprintf(stderr, "Could not open file %s for writing\n", filename);
      return(-1);
    }
    fprintf(ofs, "%4d\t%2d\t%25.16e\t%25.16e\n", nstore, dir, creal(pl), cimag(pl));
    fclose(ofs);
#if defined MPI
  }
#endif
  free(tmp_loc);
  return(0);
}
