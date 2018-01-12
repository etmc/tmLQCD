/***********************************************************************
*
* Copyright (C) 1995 Ulli Wolff, Stefan Sint
*               2001,2005 Martin Hasenbusch
*               2011,2012 Carsten Urbach
*               2013      Albert Deuzeman
*               2015,2018 Bartosz Kostrzewa
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
#ifdef TM_USE_MPI
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
#include "field_strength_types.h" 
#include "kahan_summation.h"
#include "omp_accumulator.h"
#include "tensors.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void measure_clover_field_strength_observables(const su3 ** const gf, field_strength_obs_t * const fso)
{
  // we have iG_\mu\nu = 1/4 P_T.A. [clover] where P is the projection to the
  // traceless anti-hermitian part
  // the minus sign compensates for the i^2 in the lattice definition of G_\mu\nu
  // our traceless anti-hermitian projection includes a factor of 0.5, so instead of 
  // the usual (1/8)^2 we get (1/4)^2 of the clover
  // 1/4 from the definition of the energy density <E> = 1\4 Tr[ G_\mu\nu G_\mu\nu ]
  //
  // The additional multiplication by 4 (the first factor below), originates in the fact
  // that we only accumulate over the upper triangle of G_\mu\nu below, which
  // such that iG_\mu\nu = 1/2 P_T.A. [clover upper triangle]
  
  const double energy_density_normalization = - 4 / ( 4 * 16.0 * VOLUME * g_nproc);

  // for the toplogical charge, we would naively have a normalisation of -1/(16 * 32 * pi^2)
  // but we only sum up contributions from the upper triangle of Gmunu, saving a factor of 4
  const double topo_charge_normalization = - 1 / ( 4 * 32.0 * M_PI * M_PI );
  double Eres = 0;
  double Qres = 0;
 
  // Euclidean 4D totally anti-symemtric tensor 
  epsilon4_t eps4 = new_epsilon4();

#ifdef TM_USE_OMP
  // accumulators for thread-local sums
  omp_re_acc_t Eacc = new_omp_re_acc();
  omp_re_acc_t Qacc = new_omp_re_acc();
  // this involves memory allocation, so we need corresponding frees!
  omp_re_acc_init(&Eacc, 1);
  omp_re_acc_init(&Qacc, 1);
#endif

#ifdef TM_USE_MPI
  double ALIGN mres=0;
#endif

#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif
    su3 ALIGN v1, v2, plaq;
    double ALIGN E;
    double ALIGN Q;
    
    // kahan accumulators for energy density and top. charge
    kahan_re_t E_kahan = new_kahan_re();
    kahan_re_t Q_kahan = new_kahan_re();

    // for the measurement of the top. charge density, we need to temporarily
    // store the components of Gmunu
    // for simplicity of notation, we allocate 4x4 but will only use the
    // upper triangle
    su3 Gmunu[4][4];
  
    /*  compute the clover-leaves, store them in Gmunu and compute the energy density
     *  later compute the topological charge */
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
          _su3_assign(Gmunu[k][l], plaq);
          
          // compute and accumulate the energy density at this stage
          _trace_su3_times_su3(E, plaq, plaq);
          kahan_sum_re_step(E, &E_kahan);
        }
      }
      
      // sum up the topological charge contribution now
      for( int i = 0; i < eps4.N; i++ ){
        int i1 = eps4.eps_idx[i][0];
        int i2 = eps4.eps_idx[i][1];
        int i3 = eps4.eps_idx[i][2];
        int i4 = eps4.eps_idx[i][3];

        // when Gmunu components from the lower triangle are to be used,
        // we can simply skip them and multiply our normalisation by a factor of two
        if( eps4.eps_idx[i][1] < eps4.eps_idx[i][0] ){
          continue;
        }
        if( eps4.eps_idx[i][3] < eps4.eps_idx[i][2] ){
          continue;
        }
        
        _trace_su3_times_su3( Q, 
                              Gmunu[ i1 ][ i2 ],
                              Gmunu[ i3 ][ i4 ] );

        // (Kahan) accumulate topological charge and take care of signs coming
        // the Levi-Civita
        kahan_sum_re_step(eps4.eps_val[i]*Q, &Q_kahan);
      }

    }
    
    E = kahan_sum_re_final(&E_kahan);
    Q = kahan_sum_re_final(&Q_kahan);

#ifdef TM_USE_OMP
    omp_re_acc_add(&Eacc, &E, 1);
    omp_re_acc_add(&Qacc, &Q, 1);
  } /* OpenMP parallel closing brace */

  omp_re_acc_reduce(&Eres, &Eacc);  
  omp_re_acc_reduce(&Qres, &Qacc);
  // free omp accumulator memory
  omp_re_acc_free(&Eacc);
  omp_re_acc_free(&Qacc); 

#else
  Eres = E;
  Qres = Q;
#endif

#ifdef TM_USE_MPI
  MPI_Allreduce(&Eres, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Eres = mres;
  MPI_Allreduce(&Qres, &mres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Qres = mres;
#endif
  fso->E = energy_density_normalization * Eres;
  // TODO: normalisation of top. charge
  fso->Q = topo_charge_normalization * Qres;
}
