/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2008 Carsten Urbach
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
#include <math.h>
#include <time.h>
#include "global.h"
#include "start.h"
#include "ranlxs.h"
#include "su3spinor.h"
#include "source_generation.h"
#include "invert_eo.h"
#include "solver/solver.h"
#include "geometry_eo.h"
#include "linalg/convert_eo_to_lexic.h"
#include "measurements.h"
#include "online_measurement.h"

/******************************************************
 *
 * This routine computes the correlators
 * <PP>, <PA> and <PV> (<source sink>)
 * using a stochastic time slice source
 * and only one inversion (actually A_0)
 * 
 * for <AP> we would need another inversion
 *
 *
 *
 ******************************************************/

void online_measurement(const int traj, const int id) {
  int i, j, t, tt, t0;
  double *Cpp, *Cpa, *Cp4;
  double res = 0., respa = 0., resp4 = 0.;
  double atime, etime;
  float tmp;
#ifdef MPI
  double mpi_res = 0., mpi_respa = 0., mpi_resp4 = 0.;
#endif
  FILE *ofs;
  char *filename;
  char buf[100];
  spinor phi;
  filename=buf;
  sprintf(filename,"%s%.6d", "onlinemeas." ,traj);

  /* generate random timeslice */
  if(ranlxs_init == 0) {
    rlxs_init(1, 123456);
  }
  ranlxs(&tmp, 1);
  t0 = (int)(measurement_list[id].max_source_slice*tmp);
#ifdef MPI
  MPI_Bcast(&t0, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if(g_debug_level > 1 && g_proc_id == 0) {
    printf("# timeslice set to %d (T=%d) for online measurement\n", t0, g_nproc_t*T);
    printf("# online measurements parameters: kappa = %g, mu = %g\n", g_kappa, g_mu/2./g_kappa);
  }
#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

  Cpp = (double*) calloc(g_nproc_t*T, sizeof(double));
  Cpa = (double*) calloc(g_nproc_t*T, sizeof(double));
  Cp4 = (double*) calloc(g_nproc_t*T, sizeof(double));

  source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
			      t0, 0, traj);

  invert_eo(g_spinor_field[2], g_spinor_field[3], 
	    g_spinor_field[0], g_spinor_field[1],
	    1.e-14, measurement_list[id].max_iter, CG, 1, 0, 1);

  /* now we bring it to normal format */
  /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
  convert_eo_to_lexic(g_spinor_field[DUM_MATRIX], g_spinor_field[2], g_spinor_field[3]);
  
  /* now we sum only over local space for every t */
  for(t = 0; t < T; t++) {
    j = g_ipt[t][0][0][0];
    res = 0.;
    respa = 0.;
    resp4 = 0.;
    for(i = j; i < j+LX*LY*LZ; i++) {
      res += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], g_spinor_field[DUM_MATRIX][j]);
      _gamma0(phi, g_spinor_field[DUM_MATRIX][j]);
      respa += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], phi);
      _gamma5(phi, phi);
      resp4 += _spinor_prod_im(g_spinor_field[DUM_MATRIX][j], phi);
    }

#if defined MPI
    MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    res = mpi_res;
    MPI_Reduce(&respa, &mpi_respa, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    respa = mpi_respa;
    MPI_Reduce(&resp4, &mpi_resp4, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    resp4 = mpi_resp4;
#endif
    Cpp[t+g_proc_coords[0]*T] = +res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)*2.;
    Cpa[t+g_proc_coords[0]*T] = -respa/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)*2.;
    Cp4[t+g_proc_coords[0]*T] = +resp4/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)*2.;
  }

#ifdef MPI
  /* some gymnastics needed in case of parallelisation */
  if(g_mpi_time_rank == 0) {
    MPI_Gather(&Cpp[g_proc_coords[0]*T], T, MPI_DOUBLE, Cpp, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
    MPI_Gather(&Cpa[g_proc_coords[0]*T], T, MPI_DOUBLE, Cpa, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
    MPI_Gather(&Cp4[g_proc_coords[0]*T], T, MPI_DOUBLE, Cp4, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
  }
#endif

  /* and write everything into a file */
  if(g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
    ofs = fopen(filename, "w");
    fprintf( ofs, "1  1  0  %e  %e\n", Cpp[t0], 0.);
    for(t = 1; t < g_nproc_t*T/2; t++) {
      tt = (t0+t)%(g_nproc_t*T);
      fprintf( ofs, "1  1  %d  %e  ", t, Cpp[tt]);
      tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
      fprintf( ofs, "%e\n", Cpp[tt]);
    }
    tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
    fprintf( ofs, "1  1  %d  %e  %e\n", t, Cpp[tt], 0.);

    fprintf( ofs, "2  1  0  %e  %e\n", Cpa[t0], 0.);
    for(t = 1; t < g_nproc_t*T/2; t++) {
      tt = (t0+t)%(g_nproc_t*T);
      fprintf( ofs, "2  1  %d  %e  ", t, Cpa[tt]);
      tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
      fprintf( ofs, "%e\n", Cpa[tt]);
    }
    tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
    fprintf( ofs, "2  1  %d  %e  %e\n", t, Cpa[tt], 0.);

    fprintf( ofs, "6  1  0  %e  %e\n", Cp4[t0], 0.);
    for(t = 1; t < g_nproc_t*T/2; t++) {
      tt = (t0+t)%(g_nproc_t*T);
      fprintf( ofs, "6  1  %d  %e  ", t, Cp4[tt]);
      tt = (t0+g_nproc_t*T-t)%(g_nproc_t*T);
      fprintf( ofs, "%e\n", Cp4[tt]);
    }
    tt = (t0+g_nproc_t*T/2)%(g_nproc_t*T);
    fprintf( ofs, "6  1  %d  %e  %e\n", t, Cp4[tt], 0.);
    fclose(ofs);
  }
  free(Cpp); free(Cpa); free(Cp4);
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("ONLINE: measurement done int t/s = %1.4e\n", etime - atime);
  }
  return;
}
