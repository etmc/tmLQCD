/***********************************************************************
 *
 * Copyright (C) 2011 Elena Garcia-Ramos
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
#include "global.h"
#include "start.h"
#include "su3.h"
#include "linalg_eo.h"
#include "chebyshev_polynomial_nd.h"
#include <io/eospinor.h>
#include "solver/solver.h"
#include "solver/jdher.h"
#include "solver/eigenvalues.h"
#include "X_psi.h"
#include "gamma.h"
#include "P_M_eta.h"
#include "spectral_proj.h"

double mode_n;

double mode_number(spinor * const S, double const mstarsq){

  printf("Starting mode_number calculation...\n");fflush(stdout);
  spinor **s,*s_;

  s_ = calloc(2*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(2, sizeof(spinor*));

  for(int i = 0; i < 2; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
}

  /* Computing P_M = h(X)^2 */
  h_X_sqr_eta(s[0],s[1],S,mstarsq);
  
  /* Computing the mode number  nu = (|eta>,P_M|eta>)=||h(X)^2|eta>||^2 */
  /* being |eta> the stochastic source */

  mode_n=square_norm(s[1], VOLUME, 1); 

  if(g_proc_id == 0) {
  printf("The Value of the Mode Number is %f \n", mode_n);
  }

  free(s);
  free(s_);
  return(mode_n);
}


void top_sus(spinor * const S, double const mstarsq){
  
  double mode_num, topo_sus = 0.0;
  double A = 0.0, B = 0.0, C = 0.0;
  spinor **s, *s_;

  s_ = calloc(25*VOLUMEPLUSRAND+1, sizeof(spinor));
  s  = calloc(25, sizeof(spinor*));

  for(int i = 0; i < 25; i++) {
#if (defined SSE3 || defined SSE2 || defined SSE)
    s[i] = (spinor*)(((unsigned long int)(s_)+ALIGN_BASE)&~ALIGN_BASE)+i*VOLUMEPLUSRAND;
#else
    s[i] = s_+i*VOLUMEPLUSRAND;
#endif
}

  /* s[0]=h(X)|eta>  s[2]=h(X)^2|eta>*/
  h_X_sqr_eta(s[0],s[2],S,mstarsq);

  /* s[2]=[gamma5 h(X)]|eta>*/
  gamma5(s[1],s[0], VOLUME);

  /* s[3]=[h(X) gamma5 h(X)}|eta> */
  h_X_eta(s[3], s[1], mstarsq);


  /* A = (h(X)^2|eta>,h(X)^2|eta>) */
  A=scalar_prod_r(s[2],s[2], VOLUME, 1); 

  /* B = ([h(X) gamma5 h(X)]|eta>,[h(X) gamma5 h(X)]|eta>)*/
  B=scalar_prod_r(s[3],s[3], VOLUME, 1); 
  
  /* C = ([h(X)]|eta>,[gamma5 h(X)]|eta>) */
  C=scalar_prod_r(s[0],s[1], VOLUME, 1);


  if(g_proc_id == 0) {
  printf("A = %f \n", A);
  printf("B = %f \n", B);
  printf("C = %f \n", C);
  printf("C^2 = %f \n", C*C);
  }

  free(s);
  free(s_);
}

