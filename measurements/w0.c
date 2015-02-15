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
#ifdef OMP
# include <omp.h>
#endif

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "w0.h"

#include "global.h"
#include "fatal_error.h"
#include "aligned_malloc.h"
#include "energy_density.h"
#include "expo.h"
#include "get_staples.h"
#include "gettime.h"
#include "measure_gauge_action.h"

#include "xchange/xchange_gauge.h"

void step_gradient_flow(su3** x0, su3** x1, su3** x2, su3adj** z, unsigned int type, double eps ) {
  double zfac[5] = { 0, (8.0)/(9.0), (17.0)/(36.0), (3.0)/(4.0), 1 };
  double zadjfac[3] = { 0.25, 1, 1 };
  su3** fields[4];
  su3 ALIGN w;
  su3adj z_adj;

  fields[0] = x0;
  fields[1] = x1;
  fields[2] = x2;
  fields[3] = x0;

#ifdef MPI
  xchange_gauge(x0);
#endif

  for( int f = 0; f < 3; ++f){
    for( int x = 0; x < VOLUME; ++x ){
      for( int mu = 0; mu < 4; ++mu ){
        get_staples(&w, x, mu, fields[f]);
        // implementing the Iwasaki, Symanzik or DBW2 flow from here should be a trivial extension
        if(f==0){
          _trace_lambda_mul(z[x][mu],-1,w);
        }else{
          _trace_lambda_mul(z_adj,-1,w);
          _su3adj_assign_const_times_su3adj(z_adj,zfac[2*f-1],z_adj);
          _su3adj_minus_const_times_su3adj(z_adj,zfac[2*f],z[x][mu]);
          z[x][mu] = z_adj;
        }
        _su3adj_assign_const_times_su3adj(z_adj,eps*zadjfac[f],z[x][mu]);
        exposu3(&w,&z_adj);
        _su3_times_su3(fields[f+1][x][mu],w,fields[f][x][mu]);
      }
    }
#ifdef MPI
    xchange_gauge(fields[f+1]); 
#endif
  }
  printf("%lf\n",z[2][2].d2);
}

void w0_measurement(const int traj, const int id, const int ieo) {

  double const W_target = 0.3;
  double const eps_vector[3] = {0.01,0.001,0.0001};
  
  double E[3],t[3];
  double W=0, eps;
  double t1, t2;
  
  if( g_proc_id == 0 ) {
    printf("# Doing w0 measurement.\n");
  }

  aligned_su3_field_t vt = aligned_su3_field_alloc(VOLUMEPLUSRAND);
  aligned_su3_field_t x1 = aligned_su3_field_alloc(VOLUMEPLUSRAND);
  aligned_su3_field_t x2 = aligned_su3_field_alloc(VOLUMEPLUSRAND);
  aligned_su3adj_field_t z = aligned_su3adj_field_alloc(VOLUMEPLUSRAND);

#ifdef MPI
  xchange_gauge(g_gauge_field);
#endif
  memcpy(vt.field[0],g_gauge_field[0],sizeof(su3)*4*VOLUMEPLUSRAND);

  eps = eps_vector[0];
  t[0] = 0.0;
  t[1] = 0.0;
  t[2] = 0.0;
  if(g_proc_id==0){
    for(int x = 0; x < VOLUME; ++x){
      for(int mu = 0; mu < 4; ++mu){
        if(fabs(creal(g_gauge_field[x][mu].c02)-creal(vt.field[x][mu].c02)) > 0.00001){
          printf("[%d][%d].c02: %lf %lf\n",x,mu,creal(g_gauge_field[x][mu].c02),creal(vt.field[x][mu].c02));
        }
      }
    }
    printf("plaq: %lf\n",measure_plaquette(vt.field)/(6*VOLUME*g_nproc));
  }
  t1 = gettime();
  measure_energy_density(vt.field,&E[3]);
  t2 = gettime();

  while( W_target-W > 0.000001 ) {
    if(g_proc_id==0) {
      printf("time for energy density measurement: %lf\n",t2-t1);
      printf("plaq: %lf\n",measure_plaquette(vt.field)/(6*VOLUME*g_nproc));
    }
    t[0] = t[3];
    E[0] = E[3];
    for(int step = 1; step < 3; ++step) {
      t[step] += step*eps;
      step_gradient_flow(vt.field,x1.field,x2.field,z.field,0,eps);
      t1=gettime();
      measure_energy_density(vt.field,&E[step]);
      t2=gettime();
    }
    W = t[1]*t[1]*( 2*E[1] + t[1]*((E[2]-E[0])/(2*eps)) ) ;
    
    if(g_proc_id==0){
      printf("%12.12lf %12.12lf %lf \n",E[1],W,t[1]);
    }

    if( fabs(W-W_target) < 0.01 )
      eps = eps_vector[1];
    if( fabs(W-W_target) < 0.001 )
      eps = eps_vector[2];
  }

  aligned_su3_field_free(&vt);
  aligned_su3_field_free(&x1);
  aligned_su3_field_free(&x2);
  aligned_su3adj_field_free(&z);
  
  if( g_proc_id == 0 ) {
    FILE *outfile;
    char filename[] = "w0.data";
    outfile = fopen(filename,"a");

    if( outfile == NULL ) {
      char error_message[200];
      snprintf(error_message,200,"Couldn't open %s for appending during measurement %d!",filename, id);
      fatal_error(error_message,"w0_measurement");
    }

    fprintf(outfile, "%08d %14.12lf %14.12lf %14.12lf\n",traj,E[1],W,t[1]);
    fclose(outfile);
  }

  return;
}

