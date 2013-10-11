/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch, 2012 Bartosz Kostrzewa
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
 ************************************************************************/

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#ifdef OMP
# include <omp.h>
#endif

#include <string.h>
#include <stdio.h>

#include "su3.h"
#include "geometry_eo.h"
#include "global.h"
#include "measure_oriented_plaquettes.h"
#include "fatal_error.h"


void measure_oriented_plaquettes(const su3 ** const gf, double *plaq) {
#ifdef MPI
  double ALIGN mplaq[6];
#endif

  int ix,ix1,ix2,mu1,mu2,plane;
  su3 ALIGN pr1,pr2; 
  const su3 *v,*w;
  double ALIGN pl; 
  double ALIGN ks[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double ALIGN kc[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double ALIGN tr[6],ts[6],tt[6];

  for (ix=0;ix<VOLUME;ix++){
    plane = 0;
    for (mu1=0;mu1<3;mu1++){ 
      ix1=g_iup[ix][mu1];
      for (mu2=mu1+1;mu2<4;mu2++){ 
        ix2=g_iup[ix][mu2];
        v=&gf[ix][mu1];
        w=&gf[ix1][mu2];
        _su3_times_su3(pr1,*v,*w);
        v=&gf[ix][mu2];
        w=&gf[ix2][mu1];
        _su3_times_su3(pr2,*v,*w);
        _trace_su3_times_su3d(pl,pr1,pr2);
        tr[plane]=pl+kc[plane];
        ts[plane]=tr[plane]+ks[plane];
        tt[plane]=ts[plane]-ks[plane];
        ks[plane]=ts[plane];
        kc[plane]=tr[plane]-tt[plane];
        ++plane;
      }
    }
  }

  for(int j = 0; j < 6; ++j) {
    kc[j]=(kc[j]+ks[j])/3.0;
    plaq[j] = kc[j]/(g_nproc*VOLUME);
  }

#ifdef MPI
  MPI_Allreduce(plaq, mplaq, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for(int j = 0; j < 6; j++)
    plaq[j] = mplaq[j];
#endif
  return;
}

void oriented_plaquettes_measurement(const int traj, const int id, const int ieo) {
  double plaq[6];

  if( g_proc_id == 0 ) {
    printf("# Doing oriented plaquettes measurement.\n");
  }
  measure_oriented_plaquettes((const su3** const)g_gauge_field,plaq);

  if( g_proc_id == 0 ) {
    FILE *outfile;
    char filename[] = "oriented_plaquettes.data";
    outfile = fopen(filename,"a");

    if( outfile == NULL ) {
      char error_message[200];
      snprintf(error_message,200,"Couldn't open %s for appending during measurement %d!",filename, id);
      fatal_error(error_message,"oriented_plaquettes_measurement");
    }

    fprintf(outfile, "%.8d %14.12lf %14.12lf %14.12lf %14.12lf %14.12lf %14.12lf\n",traj,plaq[0],plaq[1],plaq[2],plaq[3],plaq[4],plaq[5]);
    fclose(outfile);
  }

  return;
}

