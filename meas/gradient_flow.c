/***********************************************************************
 *
 * Copyright (C) 2013 Albert Deuzeman 
 *               2015 Bartosz Kostrzewa
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
#ifdef TM_USE_OMP
# include <omp.h>
#endif
#ifdef TM_USE_MPI
# include <mpi.h>
#endif

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "fatal_error.h"
#include "aligned_malloc.h"
#include "energy_density.h"
#include "expo.h"
#include "get_staples.h"
#include "get_rectangle_staples.h"
#include "gettime.h"
#include "measure_gauge_action.h"
#include "matrix_utils.h"
#include "xchange/xchange_gauge.h"
#include "gradient_flow.h"

void step_gradient_flow(su3 ** x0, su3 ** x1, su3 ** x2, su3 ** z, const unsigned int type, const double eps ) {
  double zfac[5] = { 1, (8.0)/(9.0), (-17.0)/(36.0), (3.0)/(4.0), -1 };
  double zepsfac[3] = { 0.25, 1, 1 };
  su3** fields[4];

  fields[0] = x0;
  fields[1] = x1;
  fields[2] = x2;
  fields[3] = x0;

#ifdef TM_USE_OMP
#pragma omp parallel
#endif
  {
 
  su3 ALIGN w,w1,w2;
  su3 ALIGN z_tmp,z_tmp1;

#ifdef TM_USE_MPI
#ifdef TM_USE_OMP
#pragma omp single
#endif
  {
  xchange_gauge(x0);
  }
#endif

  // implementation of third-order Runge-Kutta integrator following Luescher's hep-lat/1006.4518
  // this can probably be improved...

  for( int f = 0; f < 3; ++f ){
#ifdef TM_USE_OMP
#pragma omp for
#endif
    for( int x = 0; x < VOLUME; ++x ){
      for( int mu = 0; mu < 4; ++mu ){
        get_staples(&w1, x, mu, fields[f]);
        // usually we dagger the staples, but the sign convention seems to require this
        _su3_times_su3d(z_tmp,w1,fields[f][x][mu]);
        project_traceless_antiherm(&z_tmp);

        // implementing the Iwasaki, Symanzik or DBW2 flow from here should be a trivial extension
        // but it will require adding some (more) parameters and making sure that g_dbw2rand exists
        // also in the inverter if the measurement is to be carried out there
        //get_rectangle_staples_general(&w2,x,mu,fields[f]);
        //_su3_times_su3d(w1,w2,fields[f][x][mu]);

        if(f==0){
          _real_times_su3(z[x][mu],eps,z_tmp);
        }else{
          _real_times_su3(z_tmp,eps*zfac[2*f-1],z_tmp);
          _su3_refac_acc(z_tmp,zfac[2*f],z[x][mu]);
          z[x][mu] = z_tmp;
        }
        _real_times_su3(z_tmp,zepsfac[f],z[x][mu]);
        project_traceless_antiherm(&z_tmp);
        cayley_hamilton_exponent(&w,&z_tmp);
        _su3_times_su3(fields[f+1][x][mu],w,fields[f][x][mu]);
      }
    }
#ifdef TM_USE_MPI
#ifdef TM_USE_OMP
#pragma omp single
#endif
    {
    xchange_gauge(fields[f+1]); 
    }
#endif
  }
  }
}

void gradient_flow_measurement(const int traj, const int id, const int ieo) {

  double E[3],t[3], P[3];
  double W=0, eps=0.01, tsqE=0;
  double t1, t2;

  if( g_proc_id == 0 ) {
    printf("# Doing gradient flow measurement.\n");
  }
  
  FILE *outfile;
  if( g_proc_id == 0 ) {
    char filename[100];
    snprintf(filename,100,"gradflow.%06d",traj);
    outfile = fopen(filename,"w");

    if( outfile == NULL ) {
      char error_message[200];
      snprintf(error_message,200,"Couldn't open %s for writing during measurement %d!",filename, id);
      fatal_error(error_message,"gradient_flow_measurement");
    }

    fprintf(outfile, "traj t P Eplaq Esym tsqEplaq tsqEsym Wsym\n");
  }

  aligned_su3_field_t vt = aligned_su3_field_alloc(VOLUMEPLUSRAND+g_dbw2rand);
  aligned_su3_field_t x1 = aligned_su3_field_alloc(VOLUMEPLUSRAND+g_dbw2rand);
  aligned_su3_field_t x2 = aligned_su3_field_alloc(VOLUMEPLUSRAND+g_dbw2rand);
  aligned_su3_field_t z = aligned_su3_field_alloc(VOLUME);

#ifdef TM_USE_MPI
  xchange_gauge(g_gauge_field);
#endif
  memcpy(vt.field[0],g_gauge_field[0],sizeof(su3)*4*(VOLUMEPLUSRAND+g_dbw2rand));

  t[0] = E[0] = P[0] = 0.0;
  t[1] = E[1] = P[1] = 0.0;
  t[2] = E[2] = P[2] = 0.0;

  t1 = gettime();
  measure_energy_density(vt.field,&E[2]);
  P[2] = measure_plaquette(vt.field)/(6.0*VOLUME*g_nproc);
  t2 = gettime();
  if(g_proc_id==0 && g_debug_level > 2) {
    printf("time for energy density measurement: %lf\n",t2-t1);
  }

  while( t[1] < 9.99 ) {
    t[0] = t[2];
    E[0] = E[2];
    P[0] = P[2];
    for(int step = 1; step < 3; ++step) {
      t[step] = t[step-1]+eps;
      step_gradient_flow(vt.field,x1.field,x2.field,z.field,0,eps);
      measure_energy_density(vt.field,&E[step]);
      P[step] = measure_plaquette(vt.field)/(6.0*VOLUME*g_nproc);
    }
    W = t[1]*t[1]*( 2*E[1] + t[1]*((E[2]-E[0])/(2*eps)) ) ;
    tsqE = t[1]*t[1]*E[1];
    
    if(g_proc_id==0 && g_debug_level > 3){
      printf("sym(plaq)  t=%lf 1-P(t)=%1.8lf E(t)=%2.8lf(%2.8lf) t^2E=%2.8lf(%2.8lf) W(t)=%2.8lf \n",t[1],1-P[1],
        E[1],36*(1-P[1]),
        tsqE,t[1]*t[1]*36*(1-P[1]),
        W);
    }
    if(g_proc_id==0){
      fprintf(outfile,"%06d %f %2.12lf %2.12lf %2.12lf %2.12lf %2.12lf %2.12lf \n",
                      traj,t[1],P[1],
                      36*(1-P[1]),E[1],
                      t[1]*t[1]*36*(1-P[1]),tsqE,
                      W);
      fflush(outfile);
    }

  }

  aligned_su3_field_free(&vt);
  aligned_su3_field_free(&x1);
  aligned_su3_field_free(&x2);
  aligned_su3_field_free(&z);
 
  t2 = gettime();
  
  if( g_proc_id == 0 ) {
    if(g_debug_level>2){
      printf("Gradient flow measurement done in %f seconds!\n",t2-t1);
    }
    fclose(outfile);
  }

  return;
}

