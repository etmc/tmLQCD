/***********************************************************************
 *
 * Copyright (C) 2007, 2008 Jan Volkholz, Carsten Urbach
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
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "expo.h"
#include "init_stout_smear_vars.h"

su3 * gauge_field_saved;
su3 ** g_gauge_field_saved;
su3 * gauge_field_smeared;
su3 ** g_gauge_field_smeared;
su3 * C_smearing;
su3 ** g_C_smearing;
su3 * Q_smearing;
su3 ** g_Q_smearing;
su3 * Q_squared_smearing;
su3 ** g_Q_squared_smearing;
su3 * B1_smearing;
su3 ** g_B1_smearing;
su3 * B2_smearing;
su3 ** g_B2_smearing;
su3 * Gamma_smearing;
su3 ** g_Gamma_smearing;
su3 * Lambda_smearing;
su3 ** g_Lambda_smearing;

double * g_c0_smearing;
double * g_c1_smearing;

complex * g_f0_smearing;
complex * g_f1_smearing;
complex * g_f2_smearing;

complex * g_b10_smearing;
complex * g_b11_smearing;
complex * g_b12_smearing;

complex * g_b20_smearing;
complex * g_b21_smearing;
complex * g_b22_smearing;

complex * g_r10_smearing;
complex * g_r11_smearing;
complex * g_r12_smearing;

complex * g_r20_smearing;
complex * g_r21_smearing;
complex * g_r22_smearing;

su3 * stout_force_field;
su3 ** g_stout_force_field;
su3 * previous_stout_force_field;
su3 ** g_previous_stout_force_field;

/*----------------------------------------------------------------------------*/

int init_stout_smear_vars(const int V, const int stout_no_iter) 
{

  printf("Running init_stout_smear_vars\n");
  const int dim = 4 ;

  int i, k, x, mu;

  i = 0;
  k = 0;
  mu = 0;

  if (g_exposu3_no_c == 0) init_exposu3();

  /*
   *  this is the field where we store the smeared force matrices \Sigma^{(k)}_\mu(x)
   *  eqtn (44) hep-lat/0311018
   */
  gauge_field_smeared = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_gauge_field_smeared = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_gauge_field_smeared[0] = (su3*)(((unsigned long int)(gauge_field_smeared)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_gauge_field_smeared[0] = gauge_field_smeared;
#endif

  for(x = 1; x < V; x++) 
  {
    g_gauge_field_smeared[x] = g_gauge_field_smeared[x-1] + 4;
  }

  /*
   *  this is the field where we store the smeared gauge_field
   */
  gauge_field_saved = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_gauge_field_saved = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_gauge_field_saved[0] = (su3*)(((unsigned long int)(gauge_field_saved)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_gauge_field_saved[0] = gauge_field_saved;
#endif

  for(x = 1; x < V; x++) 
  {
    g_gauge_field_saved[x] = g_gauge_field_saved[x-1] + 4;
  }

  /*
   *  here we save the C matrix field from eqtn(1) in hep-lat/0311018
   */
  C_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_C_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_C_smearing[0] = (su3*)(((unsigned long int)(C_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_C_smearing[0] = C_smearing;
#endif

  for(x = 1; x < V; x++) 
  {
    g_C_smearing[x] = g_C_smearing[x-1] + 4;
  }

  /*
   *  here we save the Q matrix field from eqtn(2) in hep-lat/0311018
   */
  Q_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_Q_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_Q_smearing[0] = (su3*)(((unsigned long int)(Q_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_Q_smearing[0] = Q_smearing;
#endif

  for(x = 1; x < V; x++) 
  {
    g_Q_smearing[x] = g_Q_smearing[x-1] + 4;
  }

  /*
   *  this will hold the squared of the qbove
   */
  Q_squared_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_Q_squared_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_Q_squared_smearing[0] = (su3*)(((unsigned long int)(Q_squared_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_Q_squared_smearing[0] = Q_squared_smearing;
#endif

  for(x = 1; x < V; x++) 
  {
    g_Q_squared_smearing[x] = g_Q_squared_smearing[x-1] + 4;
  }

  /*
   *  here we save the B1  and the B2 matrix field from eqtn(69) in hep-lat/0311018
   */
  B1_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  B2_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_B1_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  g_B2_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_B1_smearing[0] = (su3*)(((unsigned long int)(B1_smearing)+ALIGN_BASE)&~ALIGN_BASE);
  g_B2_smearing[0] = (su3*)(((unsigned long int)(B2_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_B1_smearing[0] = B1_smearing;
  g_B2_smearing[0] = B2_smearing;
#endif
  for(x = 1; x < V; x++) 
  {
    g_B1_smearing[x] = g_B1_smearing[x-1] + 4;
    g_B2_smearing[x] = g_B2_smearing[x-1] + 4;
  }

  /*
   *  here we hold the Gamma matrix field from eqtn(74) in hep-lat/0311018
   */
  Gamma_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_Gamma_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_Gamma_smearing[0] = (su3*)(((unsigned long int)(Gamma_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_Gamma_smearing[0] = Gamma_smearing;
#endif

  for(x = 1; x < V; x++) 
  {
    g_Gamma_smearing[x] = g_Gamma_smearing[x-1] + 4;
  }

  /*
   *  here we save the Lambda matrix field from eqtn(73) in hep-lat/0311018
   */
  Lambda_smearing = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_Lambda_smearing = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_Lambda_smearing[0] = (su3*)(((unsigned long int)(Lambda_smearing)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_Lambda_smearing[0] = Lambda_smearing;
#endif

  for(x = 1; x < V; x++) 
  {
    g_Lambda_smearing[x] = g_Lambda_smearing[x-1] + 4;
  }


  /*
   *  these are the c_0 and c_1 fields from eqtns (14) and (15) in hep-lat/0311018
   */
  g_c0_smearing = calloc(V, sizeof(double));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  g_c1_smearing = calloc(V, sizeof(double));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  /*
   *  these are the f0, f1 and f2 fields from eqtn(29) in hep-lat/0311018
   */
  g_f0_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_f1_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_f2_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  /*
   *  these are the b10, b11f, b12, b20, b21 and b22 fields 
   *  from eqtns (57) and (58)  in hep-lat/0311018
   */
  g_b10_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_b11_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_b12_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  g_b20_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_b21_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_b22_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  /*
   *  these are the r10, r11f, r12, r20, r21 and r22 fields 
   *  from eqtns (57) and (58)  in hep-lat/0311018
   */
  g_r10_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_r11_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_r12_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  g_r20_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_r21_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_r22_smearing = calloc(V, sizeof(complex));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  /*
   *  this is the field where we store the smeared force matrices \Sigma^{(k)}_\mu(x)
   *  eqtn (44) hep-lat/0311018
   */
  stout_force_field = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_stout_force_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_stout_force_field[0] = (su3*)(((unsigned long int)(stout_force_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_stout_force_field[0] = stout_force_field;
#endif

  for(x = 1; x < V; x++) 
  {
    g_stout_force_field[x] = g_stout_force_field[x-1] + 4;
  }


  /*
   *  we need a second force field to store \Sigma'_\mu(x)
   *  eqtn (44) hep-lat/0311018
   */
  previous_stout_force_field = calloc(dim*V+1, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_previous_stout_force_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

#if (defined SSE || defined SSE2 || defined SSE3)
  g_previous_stout_force_field[0] = (su3*)(((unsigned long int)(previous_stout_force_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
  g_previous_stout_force_field[0] = previous_stout_force_field;
#endif

  for(x = 1; x < V; x++) 
  {
    g_previous_stout_force_field[x] = g_previous_stout_force_field[x-1] + 4;
  }

  /*printf("Leaving init_stout_smear_vars\n");*/
  return(0);

  /*
   *  here we save the Q matrix field from eqtn(2) in hep-lat/0311018
   */
  /*Q_smear_iterations = calloc(stout_no_iter*dim*V, sizeof(su3));
    if(errno == ENOMEM) 
    {
    return(1);
    }

    g_Q_smear_iterations = calloc(stout_no_iter, sizeof(su3**));
    if(errno == ENOMEM) 
    {
    return(1);
    }

    tmp_su3_pointer = Q_smear_iterations;
    for(i = 0; i < stout_no_iter; i++) 
    {
    g_Q_smear_iterations[i] = calloc(V, sizeof(su3*));
    if(errno == ENOMEM) 
    {
    return(1);
    }

    for(x = 0; x < V; x++)
    { 
    g_Q_smear_iterations[i][x] = tmp_su3_pointer;
    tmp_su3_pointer += dim;
    if(errno == ENOMEM) 
    {
    return(1);
    }
    }
    }*/
}

/*----------------------------------------------------------------------------*/

void free_stout_smear_vars() 
{
  free(gauge_field_saved);
  free(g_gauge_field_saved);
  free(gauge_field_smeared);
  free(g_gauge_field_smeared);
  free(C_smearing);
  free(g_C_smearing);
  free(Q_smearing);
  free(g_Q_smearing);
  free(Q_squared_smearing);
  free(g_Q_squared_smearing);
  free(B1_smearing);
  free(g_B1_smearing);
  free(B2_smearing);
  free(g_B2_smearing);
  free(Gamma_smearing);
  free(g_Gamma_smearing);
  free(Lambda_smearing);
  free(g_Lambda_smearing);
  free(g_c0_smearing);
  free(g_c1_smearing);
  free(g_f0_smearing);
  free(g_f1_smearing);
  free(g_f2_smearing);
  free(g_b10_smearing);
  free(g_b11_smearing);
  free(g_b12_smearing);
  free(g_b20_smearing);
  free(g_b21_smearing);
  free(g_b22_smearing);
  free(g_r10_smearing);
  free(g_r11_smearing);
  free(g_r12_smearing);
  free(g_r20_smearing);
  free(g_r21_smearing);
  free(g_r22_smearing);
  free(stout_force_field);
  free(g_stout_force_field);
  free(previous_stout_force_field);
  free(g_previous_stout_force_field);
}
