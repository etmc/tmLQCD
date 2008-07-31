/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

void online_measurement(const int traj, const int t0) {
  int i, j, t, tt;
  double *C;
  double res = 0.;
#ifdef MPI
  double mpi_res = 0.;
#endif
  FILE *ofs;
  char *filename;
  char buf[100];
  filename=buf;
  sprintf(filename,"%s%.6d.dat", "pp_corr." ,traj);


  C = (double*) calloc(g_nproc_t*T, sizeof(double));

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
    for(i = j; i < j+LX*LY*LZ; i++) {
      res += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], g_spinor_field[DUM_MATRIX][j]);
    }

#if defined MPI
    MPI_Allreduce(&res, &res2, 1, MPI_DOUBLE, MPI_SUM, g_mpi_time_slices);
    res = res2;
#endif
    C[t+g_proc_coords[0]*T] = res;
  }
#ifdef MPI
  if(g_mpi_time_rank == 0) {
    if(g_proc_coords[0] == 0) {
      
    }
    else {
    }
  }
#endif
  if(g_mpi_time_rank == 0 && g_proc_coords[0] == 0) {
    ofs = fopen(filename, "w");
    fprintf( ofs, "%e %e\n", C[t0], 0.);
    for(t = 1; t < g_nproc_t*T/2-1; t++) {
      tt = (t0+t)%T;
      fprintf( ofs, "%e  ", C[tt]);
      tt = (t0+g_nproc_t*T-1-t)%T;
      fprintf( ofs, "%e\n", C[tt]);
    }
    tt = (t0+g_nproc_t*T/2)%T;
    fprintf( ofs, "%e  %e\n", C[tt], 0.);
    fclose(ofs);
  }
  free(C);
  return;
}
