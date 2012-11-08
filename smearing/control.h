#ifndef GUARD_SMEARING_CONTROL_H
#define GUARD_SMEARING_CONTROL_H

#include <stdarg.h>

#include <buffers/gauge.h>
#include <buffers/adjoint.h>

#include <smearing/identity.h>
#include <smearing/ape.h>
#include <smearing/hyp.h>
#include <smearing/stout.h>
// #include <smearing/unitary.h>
// #include <smearing/hex.h>
// #include <smearing/gradient.h>
// #include <smearing/jacobi.h>


typedef enum 
{
  Identity = 0, 
  APE = 1,
  HYP = 2,
  Stout = 3,
  Unitary = 4,
  HEX = 5,
  Gradient = 6
} smearing_type;

typedef struct
{  
  smearing_type type;
  int id;

  /* Flags */
  /* FIXME These currently mirror flags in the separate control structs.
     A better solution should be found. */
  int calculate_force_terms;
  int smearing_performed;
  
  /* Results -- main output for users */
  /* Both of these are shallow copies -- essentially just pointers */
  gauge_field_t    result;
  adjoint_field_t  force_result;

  void* type_control;
} smearing_control_t;

smearing_control_t *construct_smearing_control(smearing_type type,  int calculate_force_terms, ...);
void free_smearing_control(smearing_control_t *control);

void smear(smearing_control_t *control, gauge_field_t in);
void smear_forces(smearing_control_t *control, adjoint_field_t in);

#endif
