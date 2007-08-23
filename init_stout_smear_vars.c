/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "init_stout_smear_vars.h"

su3 * gauge_field_smear_iterations;
su3 *** g_gauge_field_smear_iterations;
su3 * Q_smear_iterations;
su3 *** g_Q_smear_iterations;
su3 * C_smear_iterations;
su3 *** g_C_smear_iterations;
su3 * stout_force_field;
su3 ** g_stout_force_field;
su3 * previous_stout_force_field;
su3 ** g_previous_stout_force_field;
su3 * stout_Lambda_field;
su3 ** g_stout_Lambda_field;
su3 * stout_Gamma_field;
su3 ** g_stout_Gamma_field;

int init_stout_smear_vars(const int V, const int stout_no_iter) 
{

  printf("Running init_stout_smear_vars\n");
  const int dim = 4 ;

  int i, k, x, mu;

  su3* tmp_su3_pointer;
  i = 0;
  k = 0;
  mu = 0;

  /*
   *  this is the field where we store the intermediate smeared gauge configurations U^{(i)}
   *  eqtn (4) hep-lat/0311018
   */
  gauge_field_smear_iterations = calloc((stout_no_iter+1)*dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  g_gauge_field_smear_iterations = calloc(stout_no_iter+1, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = gauge_field_smear_iterations;
  for(i = 0; i < stout_no_iter+1; i++) 
  {
    g_gauge_field_smear_iterations[i] = calloc(V, sizeof(su3*));
    if(errno == ENOMEM) 
    {
      return(1);
    }

    for(x = 0; x < V; x++)
    { 
      g_gauge_field_smear_iterations[i][x] = tmp_su3_pointer;
      tmp_su3_pointer += dim;
      /*g_gauge_field_smear_iterations[i][x] = &(gauge_field_smear_iterations[i*V*dim+x*dim]);*/
      if(errno == ENOMEM) 
      {
        return(1);
      }
    }
  }

  /*
   *  here we save the C matrix field from eqtn(1) in hep-lat/0311018
   */
  C_smear_iterations = calloc(stout_no_iter*dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  
  g_C_smear_iterations = calloc(stout_no_iter, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = C_smear_iterations;
  for(i = 0; i < stout_no_iter; i++) 
  {
    g_C_smear_iterations[i] = calloc(V, sizeof(su3*));
    if(errno == ENOMEM) 
    {
      return(1);
    }

    for(x = 0; x < V; x++)
    { 
      g_C_smear_iterations[i][x] = tmp_su3_pointer;
      tmp_su3_pointer += dim;
      if(errno == ENOMEM) 
      {
        return(1);
      }
    }
  }

  /*
   *  here we save the Q matrix field from eqtn(2) in hep-lat/0311018
   */
  Q_smear_iterations = calloc(stout_no_iter*dim*V, sizeof(su3));
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
  }

  /*
   *  this is the field where we store the smeared force matrices \Sigma^{(k)}_\mu(x)
   *  eqtn (44) hep-lat/0311018
   */
  stout_force_field = calloc(dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = stout_force_field;
  g_stout_force_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  for(x = 0; x < V; x++) 
  {
    g_stout_force_field[x] = tmp_su3_pointer;
    tmp_su3_pointer += dim;
    if(errno == ENOMEM) 
    {
      return(1);
    }
  }

  /*
   *  we need a second force field to store \Sigma'_\mu(x)
   *  eqtn (44) hep-lat/0311018
   */
  previous_stout_force_field = calloc(dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = previous_stout_force_field;
  g_previous_stout_force_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  for(x = 0; x < V; x++) 
  {
    g_previous_stout_force_field[x] = tmp_su3_pointer;
    tmp_su3_pointer += dim;
    if(errno == ENOMEM) 
    {
      return(1);
    }
  }

  /*
   *  this is the field where we store the \Lambda^{(k)}_\mu(x) 
   *  eqtn (73) hep-lat/0311018
   */
  stout_Lambda_field = calloc(dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = stout_Lambda_field;
  g_stout_Lambda_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  for(x = 0; x < V; x++) 
  {
    g_stout_Lambda_field[x] = tmp_su3_pointer;
    tmp_su3_pointer += dim;
    if(errno == ENOMEM) 
    {
      return(1);
    }
  }

  /*
   *  this is the field where we store the \Gamma^{(k)}_\mu(x) 
   *  eqtn (73) hep-lat/0311018
   */
  stout_Gamma_field = calloc(dim*V, sizeof(su3));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  tmp_su3_pointer = stout_Gamma_field;
  g_stout_Gamma_field = calloc(V, sizeof(su3**));
  if(errno == ENOMEM) 
  {
    return(1);
  }

  for(x = 0; x < V; x++) 
  {
    g_stout_Gamma_field[x] = tmp_su3_pointer;
    tmp_su3_pointer += dim;
    if(errno == ENOMEM) 
    {
      return(1);
    }
  }

  /*
   *  here we allocate the testspinors necessary to get the 
   *  explicit force \Sigma
   */
  /*g_test_spinor_field_left = calloc(V/2, sizeof(spinor));
  if(errno == ENOMEM) 
  {
    return(1);
  }
  
  g_test_spinor_field_right = calloc(V/2, sizeof(spinor));
  if(errno == ENOMEM) 
  {
    return(1);
  }*/
  
  /*#if (defined SSE || defined SSE2 || defined SSE3)
    g_gauge_field[0] = (su3*)(((unsigned long int)(gauge_field)+ALIGN_BASE)&~ALIGN_BASE);
#else
g_gauge_field[0] = gauge_field;
#endif*/
  /*for(i = 1; i < V; i++){
    g_gauge_field[i] = g_gauge_field[i-1]+4;
    }*/

  /*#  if defined _USE_HALFSPINOR
    if(back == 1) {*/
  /*
     g_gauge_field_copy[ieo][PM][sites/2][mu]
   */
  /*    g_gauge_field_copy = (su3***)calloc(2, sizeof(su3**));
        g_gauge_field_copy[0] = (su3**)calloc(VOLUME, sizeof(su3*));
        g_gauge_field_copy[1] = g_gauge_field_copy[0] + (VOLUME)/2;
        if(errno == ENOMEM) {
        return(2);
        }
        gauge_field_copy = (su3*)calloc(4*(VOLUME)+1, sizeof(su3));
        if(errno == ENOMEM) {
        return(2);
        }
#    if (defined SSE || defined SSE2 || defined SSE3)
g_gauge_field_copy[0][0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#    else
g_gauge_field_copy[0][0] = gauge_field_copy;
#    endif
for(i = 1; i < (VOLUME)/2; i++) {
g_gauge_field_copy[0][i] = g_gauge_field_copy[0][i-1]+4;
}
g_gauge_field_copy[1][0] = g_gauge_field_copy[0][0] + 2*VOLUME; 
for(i = 1; i < (VOLUME)/2; i++) {
g_gauge_field_copy[1][i] = g_gauge_field_copy[1][i-1]+4;
}
}
#  else
if(back == 1) {
g_gauge_field_copy = calloc((VOLUME+RAND), sizeof(su3*));
if(errno == ENOMEM) {
return(2);
}
gauge_field_copy = calloc(8*(VOLUME+RAND)+1, sizeof(su3));
if(errno == ENOMEM) {
return(2);
}
#  if (defined SSE || defined SSE2 || defined SSE3)
g_gauge_field_copy[0] = (su3*)(((unsigned long int)(gauge_field_copy)+ALIGN_BASE)&~ALIGN_BASE);
#  else
g_gauge_field_copy[0] = gauge_field_copy;
#  endif
for(i = 1; i < (VOLUME+RAND); i++) {
g_gauge_field_copy[i] = g_gauge_field_copy[i-1]+8;
}
}
#  endif*/

  /*printf("Leaving init_stout_smear_vars\n");*/
  return(0);
  }

void free_stout_smear_vars() 
{
  free(gauge_field_smear_iterations);
  free(g_gauge_field_smear_iterations);
  free(C_smear_iterations);
  free(g_C_smear_iterations);
  free(Q_smear_iterations);
  free(g_Q_smear_iterations);
  free(stout_force_field);
  free(g_stout_force_field);
  free(stout_Lambda_field);
  free(g_stout_Lambda_field);
  free(stout_Gamma_field);
  free(g_stout_Gamma_field);
}
