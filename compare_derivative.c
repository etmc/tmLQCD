/***********************************************************************
 *
 * Copyright (C) 2024 Bartosz Kostrzewa
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
# include<tmlqcd_config.h>
#endif
#ifdef TM_USE_OMP 
# include <omp.h>
#endif
#include <stdio.h>
#include "global.h"
#include "monomial/monomial.h"

/* this function compares two derivatives calculated by an external library and tmLQCD */
void compare_derivative(monomial *mnl, su3adj **ext_lib, su3adj **native, 
    const double threshold, const char * name){
  int n_diff = 0;

  for(int ix = 0; ix < VOLUME; ix++){
    for(int mu=0; mu<4; mu++){
      double *ext=&(ext_lib[ix][mu].d1);
      double *nat=&(native[ix][mu].d1);
      for(int j=0; j<8; ++j){
        double diff=ext[j]-nat[j];
        if (sqrt(diff*diff) > threshold || isnan( ext[j] ) || isinf(ext[j]) ){
            n_diff++;
            printf("derivative at (t,x,y,z,mu,j) %d,%d,%d,%d,%d,%d,"
                   " ext: %-14e, native: %-14e ratio: %-14g diff %-14g  on proc_id %d\n", 
                  g_coord[ix][0], g_coord[ix][1], g_coord[ix][2], g_coord[ix][3], mu, j,
                  ext[j], nat[j], ext[j]/nat[j], ext[j]-nat[j], g_proc_id);
        }
      }
    }
  }
  if(n_diff > 0){
    printf("%s: the deviation between tmLQCD and the external library "
           "exceeds the threshold %.1e in %d case(s) for parameters: c0=%e c1=%e g_beta=%e on proc_id: %d\n",
           name,
           threshold,
           n_diff,
           mnl->c0,
           mnl->c1,
           mnl->beta,
           g_proc_id);

    if(g_strict_residual_check) fatal_error("Difference between external library and tmLQCD-native function!", 
                                            name);
  }
}

