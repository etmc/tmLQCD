/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
 *               2009 Florian Burger
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
#include "ranlxd.h"
#include "ranlxs.h"
#include "su3spinor.h"
#include "source_generation.h"
#include "invert_eo.h"
#include "solver/solver.h"
#include "solver/solver_params.h"
#include "geometry_eo.h"
#include "linalg/convert_eo_to_lexic.h"
#include "measurements.h"
#include "pion_norm.h"
#include "gettime.h"

void pion_norm(const int traj, const int id, const int ieo) {
  int i, j, z, zz, z0;
  double *Cpp;
  double res = 0.;
  double pionnorm;
  double atime, etime;
  float tmp;
  solver_params_t tmp_solver_params;
#ifdef MPI
  double mpi_res = 0.;
#endif
  FILE *ofs, *ofs2;
  char *filename, *filename2, *sourcefilename;
  char buf[100];
  char buf2[100];
  char buf3[100];
  filename=buf;
  filename2=buf2;
  sourcefilename=buf3;
  sprintf(filename,"pionnormcorrelator_finiteT.%.6d",traj);
  sprintf(filename2,"%s", "pion_norm.data");

  /* generate random source point */
  if(ranlxs_init == 0) {
    rlxs_init(1, 123456);
  }
  ranlxs(&tmp, 1);
  z0 = (int)(measurement_list[id].max_source_slice*tmp);
#ifdef MPI
  MPI_Bcast(&z0, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  atime = gettime();

  Cpp = (double*) calloc(g_nproc_z*LZ, sizeof(double));

  printf("Doing finite Temperature online measurement\n");
  
  /* stochastic source in z-slice */
  source_generation_pion_zdir(g_spinor_field[0], g_spinor_field[1], 
                            z0, 0, traj);
  

  /* invert on the stochastic source */
  invert_eo(g_spinor_field[2], g_spinor_field[3], 
            g_spinor_field[0], g_spinor_field[1],
            1.e-14, measurement_list[id].max_iter, CG, 1, 0, ieo, 0, NULL,tmp_solver_params, -1);

  /* now we bring it to normal format */
  /* here we use implicitly DUM_MATRIX and DUM_MATRIX+1 */
  convert_eo_to_lexic(g_spinor_field[DUM_MATRIX], g_spinor_field[2], g_spinor_field[3]);
  
  /* now we sums only over local space for every z */
  for(z = 0; z < LZ; z++) {
    res = 0.;
    /* sum here over all points in one z-slice 
       we have to look up g_ipt*/

    j = g_ipt[0][0][0][z];
    for(i = 0; i < T*LX*LY ; i++) {
           res += _spinor_prod_re(g_spinor_field[DUM_MATRIX][j], g_spinor_field[DUM_MATRIX][j]);
           j += LZ; /* jump LZ sites in array, z ist fastest index */
    }


    
#if defined MPI
    MPI_Reduce(&res, &mpi_res, 1, MPI_DOUBLE, MPI_SUM, 0, g_mpi_z_slices);
    res = mpi_res;
#endif
    Cpp[z+g_proc_coords[3]*LZ] = +res/(g_nproc_x*LX)/(g_nproc_y*LY)/(g_nproc_t*T)*2.;
  }

#ifdef MPI
  /* some gymnastics needed in case of parallelisation */
  if(g_mpi_z_rank == 0) {
    MPI_Gather(&Cpp[g_proc_coords[3]*LZ], LZ, MPI_DOUBLE, Cpp, LZ, MPI_DOUBLE, 0, g_mpi_ST_slices);
  }
#endif


  /* and write everything into a file */
  if(g_mpi_z_rank == 0 && g_proc_coords[3] == 0) {
    ofs = fopen(filename, "w");
    fprintf( ofs, "1  1  0  %e  %e\n", Cpp[z0], 0.);
    for(z = 1; z < g_nproc_z*LZ/2; z++) {
      zz = (z0+z)%(g_nproc_z*LZ);
      fprintf( ofs, "1  1  %d  %e  ", z, Cpp[zz]);
      zz = (z0+g_nproc_z*LZ-z)%(g_nproc_z*LZ);
      fprintf( ofs, "%e\n", Cpp[zz]);
    }
    zz = (z0+g_nproc_z*LZ/2)%(g_nproc_z*LZ);
    fprintf( ofs, "1  1  %d  %e  %e\n", z, Cpp[zz], 0.);
    fclose(ofs);
    
    /* sum over all Cpp to get pionnorm*/
    ofs2 = fopen(filename2, "a");
    pionnorm = 0.;
    for(z=0; z<g_nproc_z*LZ; z++){
      pionnorm += Cpp[z];
    }
    /* normalize */
    pionnorm = pionnorm/(g_nproc_z*LZ); 
    fprintf(ofs2,"%d\t %.16e\n",traj,pionnorm);
    fclose(ofs2);
  }
  
  free(Cpp);
  etime = gettime();
  if(g_proc_id == 0 && g_debug_level > 0) {
    printf("PIONNORM : measurement done int t/s = %1.4e\n", etime - atime);
  }
  return;
}
/*end  Florian Burger 4.11.2009 */

