#pragma once

#include <smearing/stout.h>
#include <smearing/hyp.h>

typedef su3 su3_two_index[12];
typedef stout_notes_t stout_notes_two_index[12];

typedef su3 su3_two_index[24];
typedef stout_notes_t stout_notes_two_index[24];

typedef struct
{
  /* Flags */
  int             calculate_force_terms;
  int             smearing_performed;
  
  /* Parameters */
  double          rho[3];
  unsigned int    iterations;
   
  /* Results -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */
  adjoint_field_t  force_result;
  
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t      *U;     /* The sequence of iterations gauge fields */
   
  su3_outer         **V_stage_0;
  stout_notes_outer **trace_stage_0;
  
  su3_outer         **V_stage_1;
  stout_notes_outer **trace_stage_1;
  
  stout_notes_tuple **trace_stage_2; /* Intermediate results to avoid double calculations */
  
  /* Final results -- the first is a shallow copy */
  gauge_field_t    result;
  adjoint_field_t  force_result;
} hex_control;

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const rho_1, double const rho_2, double const rho_3);
void free_hex_control(hex_control *control);

void hex_smear(hex_control *control, gauge_field_t in);
void hex_smear_forces(hex_control *control, adjoint_field_t in);
