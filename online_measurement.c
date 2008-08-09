/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "start.h"
#include "ranlxd.h"
#include "su3spinor.h"
#include "source_generation.h"
#include "invert_eo.h"
#include "solver/solver.h"
#include "geometry_eo.h"
#include "linalg/convert_eo_to_lexic.h"
#include "online_measurement.h"

/******************************************************
 *
 * This routine computes the two correlator
 * <PP> and <PA>  (<source sink>)
 * using a stochastic time slice source
 * and only one inversion (actually A_0)
 * 
 * for <AP> we would need another inversion
 *
 *
 *
 ******************************************************/

void online_measurement(const int traj, const int t0) {
  int i, j, t, tt;
  double *Cpp, *Cpa;
  double res = 0., respa = 0.;
  double atime, etime;
#ifdef MPI
  double mpi_res = 0., mpi_respa = 0.;
#endif
  FILE *ofs;
  char *filename;
  char buf[100];
  spinor phi;
  filename=buf;
  sprintf(filename,"%s%.6d", "onlinemeas." ,traj);

#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = (double)clock()/(double)(CLOCKS_PER_SEC);
#endif

  Cpp = (double*) calloc(g_nproc_t*T, sizeof(double));
  Cpa = (double*) calloc(g_nproc_t*T, sizeof(double));

  source_generation_pion_only(g_spinor_field[0], g_spinor_field[1], 
			      t0, 0, traj);
  
  invert_eo(g_spinor_field[2], g_spinor_field[3], 
	    g_spinor_field[0], g_spinor_field[1],
	    1.e-14, 10000, CG, 1, 0, 1);

  /* now we bring it to normal format */
  /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
  convert_eo_to_lexic(g_spinor_field[DUM_MATRIX], g_spinor_field[2], g_spinor_field[3]);
  
  /* now we sums only over local space for every t */
  for(t = 0; t < T; t++) {
    j = g_ipt[t][0][0][0];
    res = 0.;
    respa = 0.;
    for(i = j; i < j+LX*LY*LZ; i++) {
      res += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], g_spinor_field[DUM_MATRIX][j]);
      _gamma0(phi, g_spinor_field[DUM_MATRIX][j]);
      respa += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], phi);
    }

#if defined MPI
    MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    res = mpi_res;
    MPI_Reduce(&respa, &mpi_respa, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_time_slices);
    respa = mpi_respa;
#endif
    Cpp[t+g_proc_coords[0]*T] = res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/g_kappa/g_kappa/2.;
    Cpa[t+g_proc_coords[0]*T] = respa/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_z*LZ)/g_kappa/g_kappa/2.;
  }
#ifdef MPI

  /* some gymnastics needed in case of parallelisation */
  if(g_mpi_time_rank == 0) {
    MPI_Gather(&Cpp[g_proc_coords[0]*T], T, MPI_DOUBLE, Cpp, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
    MPI_Gather(&Cpa[g_proc_coords[0]*T], T, MPI_DOUBLE, Cpa, T, MPI_DOUBLE, 0, g_mpi_SV_slices);
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
    fclose(ofs);
  }
  free(Cpp); free(Cpa);
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
