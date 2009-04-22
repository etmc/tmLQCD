/*******************************************
*
* FILE: sf_observables.c
*
* Author: Jenifer Gonzalez Lopez
*
********************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "geometry_eo.h"
#include "sf_calc_action.h"
#include "sf_observables.h"

void sf_observables() {

  double plaquette_energy;
  double rectangle_energy;
  double wilson_action;
  double wilson_action_sepbound;
  double iwasaki_action;
  double partial_iwa;
  double partial_iwasaki_action;
  //su3 ** sf_background_gauge_field = NULL;

  /*** ASIGN MEMORY to the back gauge field ***/


  /*** INITIALISE the background gauge field ***/



  /* sf b.c. abelian field and standard sf weight factors included (only plaquette here) */
  plaquette_energy = measure_plaquette_sf_weights(g_Tbsf);
  wilson_action = measure_wilson_action_sf(g_Tbsf, g_beta);
  wilson_action_sepbound = measure_wilson_action_sf_separate_boundary(g_Tbsf, g_beta);
  if(g_proc_id==0){
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and standard sf weight factors included (only plaquette): \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Wilson action value sep bound is %e\n", wilson_action_sepbound); fflush(stdout); 
  }
  /* sf b.c. abelian field and weight factors for O(a)-improvement included (only plaquette here) */
  plaquette_energy = measure_plaquette_sf_weights_improvement(g_Tbsf, g_Cs, g_Ct) ;
  wilson_action = measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct);
  wilson_action_sepbound = measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct);
  if(g_proc_id==0){
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and weight factors for O(a)-improvement included (only plaquette): \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Wilson action value sep bound is %e\n", wilson_action_sepbound); fflush(stdout);
  }    
  /* sf b.c. abelian field and weight factors for O(a)-improvement included (plaquette and rectangle) */
  plaquette_energy = measure_plaquette_sf_iwasaki(g_Tbsf, g_Cs, g_Ct, g_rgi_C0) ;
  rectangle_energy = measure_rectangle_sf_iwasaki(g_Tbsf, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  if(g_proc_id==0){    
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and weight factors for O(a)-improvement included (Iwasaki = plaquette and rectangle): \n");
    fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);
  }  

  
#if 0
  /*** CHECKS: here we calculate "S[V], Gamma[V], S'[V] and Gamma'[V]" in two ways and both should agree ***/
  printf("\n"); fflush(stdout);
  printf("CHECKS: \n"); fflush(stdout);

  /* (0): here we have not yet assigned: U = V forall x_0 ==> it should not agree with the (1) and (2) */
  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action = partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_eta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf("\n"); fflush(stdout);
  printf(" Before assigning U=V but SF b.c. \n"); fflush(stdout);
  printf("S[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout); 
  printf("\n"); fflush(stdout);
  
  /* (1): identifying the gauge fields "g_gauge_fields = V" and then calculating the plaquette as usually */
  induced_lattice_background(g_gauge_field, g_Tbsf, g_eta);
  
  wilson_action = measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct);
  wilson_action_sepbound = measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct);
  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action = partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf(" Assigning U=V with the functions defined for that and then calculating S[V] from the same functions to calculate the actions as in previous cases \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", wilson_action_sepbound); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", wilson_action ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("G[V] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("S'[V] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("G'[V] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);


  /* obtain normalization factor by calculation Wilson action for U=1 in all the lattice
   and substract it to the previous result for the action.
  Therefore, it should agree with the result obtained from the analytical expression implemented below */
  set_all_links_to_one_with_dirichlet(g_Tbsf);

  iwasaki_action -= measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action -= partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf("\n"); fflush(stdout);
  printf(" Previous case but substracting the normalization factor to the action: \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("Norm - S_sf_iwasaki_notsepb[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("Norm - G[V] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("Norm' - S'[V] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("Norm' - G'[V] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout);
  printf("\n"); fflush(stdout);


  /* (2): directly from the analytical expression which has been implemente in: */
  printf("\n"); fflush(stdout);
  printf(" Assigning U=V: but directly using the analytical expression of the action S[V] \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("S[V]_analy = %e \n", lattice_background_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("G[V]_analy = %e \n", lattice_lo_effective_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("S'[V]_analy = %e \n", partial_lattice_background_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("G'[V]_analy = %e \n", partial_lattice_lo_effective_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("\n"); fflush(stdout);


  /* obtaine normalization factor by calculation Wilson action for U=1 in all the lattice */
  set_all_links_to_one_with_dirichlet(g_Tbsf);

  printf("\n"); fflush(stdout);
  printf(" Setting U=Id and Dirichlet at x0= 0, t \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  /* The next three prints give me the same result, from 3 different functions.
   The first two functions were cross-checked with Dru ==> they should be right.
   Hoever, the result here obtained still differs to what we obtain by doing the
   differenct between our result (for U=V) and the analytical expression */
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct)); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct) ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts));fflush(stdout); 
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);


  /* obtaine normalization factor by calculation Wilson action for U=1 in all the lattice */
  set_all_links_to_one();

  printf("\n"); fflush(stdout);
  printf(" Setting U=Id \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  /* For the first case below, pbc, I've gotten the number I expected: "(Nc*12*L^4)/g02".
   Thus, since the function "measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1))" was crosschecked bf with Dru it should be right.
  It somehow tells me that also the function which assigns the gauge fields to one "set_all_links_to_one()" should be right.*/
  printf("S_pbc[U,W',W] = %e \n", measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1)); fflush(stdout);
  /* The next three prints give me the same result, from 3 different functions.
   The first two functions were cross-checked with Dru ==> they should be right.
   Hoever, the result here obtained still differs to what we obtain by doing the
   differenct between our result (for U=V) and the analytical expression */
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct)); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct) ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts));fflush(stdout); 
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);

#endif

  /*** FREE MEMORY of the back gauge field ***/
  //free(sf_background_gauge_field);
}
