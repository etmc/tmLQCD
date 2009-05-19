/*******************************************
*
* FILE: sf_calc_action.h
*
* Author: Jenifer Gonzalez Lopez
*
********************************************/
#ifndef _SF_CALC_ACTION_H
#define _SF_CALC_ACTION_H

void dirichlet_boundary_conditions(int t);
void dirichlet_boundary_conditions_spatial_links_to_one(int t);
void set_all_links_to_one();
void set_all_links_to_one_with_dirichlet(int t);
void print_su3_matrix (su3 u);
void sf_boundary_conditions_spatially_constant_abelian_field(int t, double eta);

/*** MEASUREMENTS ***/

double measure_plaquette();
double measure_rectangle();

double measure_plaquette_sf_weights(int t);
double measure_plaquette_sf_weights_improvement(int t, double cs, double ct);
double measure_plaquette_sf_weights_bulk(int t);
double measure_plaquette_sf_weights_boundary_0 ();
double measure_plaquette_sf_weights_boundary_t (int t);

double measure_plaquette_sf_weights_improved_bulk(int t);
double measure_plaquette_sf_weights_improved_boundary_0 (double cs, double ct);
double measure_plaquette_sf_weights_improved_boundary_t (int t, double cs);
double measure_plaquette_sf_weights_improved_boundary_t_minus_1 (int t, double ct);

double measure_plaquette_sf_iwasaki(int t, double cs, double ct, double c0);
double measure_rectangle_sf_iwasaki(int t, double c1, double c1_ss, double c1_tss, double c1_tts);

/***** ACTIONS *****/

/*** for PBC ***/
double measure_wilson_action(double beta);
double measure_iwasaki_action(double beta, double c0, double c1);

/*** SF boundary conditions ***/
double measure_wilson_action_sf(int t, double beta);
double measure_wilson_action_sf_weights_improvement(int t, double beta, double cs, double ct);
double measure_wilson_action_sf_separate_boundary(int t, double beta);
double measure_wilson_action_sf_weights_improvement_separate_boundary(int t, double beta, double cs, double ct);
double measure_iwasaki_action_sf(int t, double beta, double cs, double ct, double c0, double c1, double c1_ss, double c1_tss, double c1_tts);

/*** FUNCTIONS NEEDED FOR THE BACKGROUND FIELD ACTION and
     BACKGROUND FIELD ACTION AND DERIVATIVE WITH RESPECT TO ETA ***/

/** PLAQUETTE (only) **/
void induced_continuum_background(su3 **b, int t, double eta);
void induced_lattice_background(su3 **v, int t, double eta);
double lattice_background_plaquette_action_sf(int t, double beta, double ct, double eta);
double lattice_lo_effective_plaquette_action_sf(int t, double beta, double ct, double eta);
double partial_lattice_background_plaquette_action_sf(int t, double beta, double ct, double eta);
double partial_lattice_lo_effective_plaquette_action_sf(int t, double beta, double ct, double eta);
double partial_lattice_lo_effective_plaquette_action_sf_k(int t, double beta, double ct, double eta);
/** IWASAKI **/
double partial_lattice_lo_effective_iwasaki_action_sf_k(int t, double beta, double c0, double c1, double eta);


/*** DEFINITION OF THE RUNNING COUPLING ***/
double partial_plaquette_sf_respect_to_eta(int t, double ct);
double partial_rectangle_sf_respect_to_eta(int t, double c1_tss, double c1_tts);
double partial_wilson_action_sf_respect_to_eta(int t, double beta, double cs, double ct);
double partial_iwasaki_action_sf_respect_to_eta(int t, double beta, double cs, double ct, double c0, double c1, double c1_ss, double c1_tss, double c1_tts);

#endif
