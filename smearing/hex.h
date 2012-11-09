#pragma once

#include <buffers/gauge.h>

typedef struct
{
  su3    Q;
  su3    expiQ;
  su3    B1;
  su3    B2;
  double f1;
  double f2;
} hex_notes_t;

typedef hex_notes_t hex_notes[4];
typedef su3 su3_outer[12];
typedef hex_notes_t hex_notes_outer[12];

typedef struct
{
  double alpha[3];
  unsigned int    iterations;
  
  /* Flags */
  int             calculate_force_terms;
  int             smearing_performed;
  
  /* Results -- main output for users */
  gauge_field_t    result; /* For direct access to the result, shallow copy... */
  adjoint_field_t  force_result;
  
  /* Intermediate results, stored to enhance locality of the analysis */
  gauge_field_t    *U;     /* The sequence of iterations gauge fields */
  hex_notes **trace; /* Intermediate results to avoid double calculations */
  
  su3_outer       **U_outer;     /* The sequence of iterations gauge fields */
  hex_notes_outer **trace_outer; /* Intermediate results to avoid double calculations */
  
} hex_control;

hex_control *construct_hex_control(int calculate_force_terms, unsigned int iterations, double const alpha_1, double const alpha_2, double const alpha_3);
void free_hex_control(hex_control *control);

void hex_smear(hex_control *control, gauge_field_t in);
void hex_smear_forces(hex_control *control, adjoint_field_t in);