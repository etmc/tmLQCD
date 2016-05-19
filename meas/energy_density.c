/***********************************************************************
*
* Copyright (C) 1995 Ulli Wolff, Stefan Sint
*               2001,2005 Martin Hasenbusch
*               2011,2012 Carsten Urbach
*               2013      Albert Deuzeman
*               2015      Bartosz Kostrzewa
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
***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#ifdef TM_USE_OMP
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "su3adj.h"
#include "matrix_utils.h"

void measure_energy_density(const su3 ** const gf, double *ret)
{
  // we have iG_\mu\nu = 1/4 P_T.A. [clover] where P is the projection to the
  // traceless anti-hermitian part
  // the minus sign compensates for the i^2 in the lattice definition of G_\mu\nu
  // our traceless anti-hermitian projection includes a factor of 0.5, so instead of 
  // the usual (1/8)^2 we get (1/4)^2 of the clover
  // 1/4 from the definition of the energy density <E> = 1\4 (G_\mu\nu)^2
  // The factor of 4 makes the result agree (at large t and keeping in mind discretization errors)
  //  with the plaquette definition and with papers... I don't understand where it comes from...
  double normalization = - 4 / ( 4 * 16.0 * VOLUME * g_nproc);
  double res = 0;
#ifdef MPI
  double ALIGN mres=0;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
    su3 ALIGN v1, v2, plaq;
    double ALIGN ac,tr,ts,tt,kc=0,ks=0;
    su3 ALIGN trace;
  
    /*  compute the clover-leave */
  /*  l  __   __
        |  | |  |
        |__| |__|
        __   __
        |  | |  |
        |__| |__| k  */
  
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for(int x = 0; x < VOLUME; x++)
    {
      for(int k = 0; k < 4; k++)
      {
        for(int l = k+1; l < 4; l++)
        {
          int xpk = g_iup[x][k];
          int xpl = g_iup[x][l];
          int xmk = g_idn[x][k];
          int xml = g_idn[x][l];
          int xpkml = g_idn[xpk][l];
          int xplmk = g_idn[xpl][k];
          int xmkml = g_idn[xml][k];
          const su3 *w1 = &gf[x][k];
          const su3 *w2 = &gf[xpk][l];
          const su3 *w3 = &gf[xpl][k];
          const su3 *w4 = &gf[x][l];
          _su3_times_su3(v1, *w1, *w2);
          _su3_times_su3(v2, *w4, *w3);
          _su3_times_su3d(plaq, v1, v2);
          w1 = &gf[x][l];
          w2 = &gf[xplmk][k];
          w3 = &gf[xmk][l];
          w4 = &gf[xmk][k];
          _su3_times_su3d(v1, *w1, *w2);
          _su3d_times_su3(v2, *w3, *w4);
          _su3_times_su3_acc(plaq, v1, v2);
          w1 = &gf[xmk][k];
          w2 = &gf[xmkml][l];
          w3 = &gf[xmkml][k];
          w4 = &gf[xml][l];
          _su3_times_su3(v1, *w2, *w1);
          _su3_times_su3(v2, *w3, *w4);
          _su3d_times_su3_acc(plaq, v1, v2);
          w1 = &gf[xml][l];
          w2 = &gf[xml][k];
          w3 = &gf[xpkml][l];
          w4 = &gf[x][k];
          _su3d_times_su3(v1, *w1, *w2);
          _su3_times_su3d(v2, *w3, *w4);
          _su3_times_su3_acc(plaq, v1, v2);
          project_traceless_antiherm(&plaq);  
          _trace_su3_times_su3(ac, plaq, plaq); // This should actually be the energy density already...
          
          // Kahan summation for each thread
          tr=ac+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
        }
      }   
    }
    kc=kc+ks;
#ifdef TM_USE_OMP
    int thread_num = omp_get_thread_num();
    g_omp_acc_re[thread_num] = kc;
  } /* OpenMP parallel closing brace */

  for(int i=0; i < omp_num_threads; ++i) {
    res += g_omp_acc_re[i];
  }
#else
  res = kc;
#endif
#ifdef MPI
  MPI_Allreduce(&res, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  res = mres;
#endif
  *ret = normalization * res;
}
