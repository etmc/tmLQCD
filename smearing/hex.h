#pragma once

#include <smearing/stout.h>
#include <smearing/hyp.h>

typedef stout_notes_t stout_notes_two_index[12];
typedef stout_notes_t stout_notes_three_index[24];

typedef struct
{
  /* Flags */
  int             calculate_force_terms;
  
  /* Parameters */
  double          coef[3];
  unsigned int    iterations;
  
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t *U;       /* The sequence of iterations gauge fields */
  
  su3_two_index           **V_stage_1;
  stout_notes_three_index **trace_stage_1;
  
  su3_two_index           **V_stage_2;
  stout_notes_two_index   **trace_stage_2;

  stout_notes_tuple       **trace_stage_3; /* Intermediate results to avoid double calculations */
  
  /* Final results -- the first is a shallow copy */
  gauge_field_t    result;
  adjoint_field_t  force_result;
} hex_control;

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const coef_0, double const coef_1, double const coef_2);
void free_hex_control(hex_control *control);

void hex_smear(hex_control *control, gauge_field_t in);
void hex_smear_forces(hex_control *control, adjoint_field_t in);
