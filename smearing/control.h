#ifndef GUARD_SMEARING_CONTROL_H
#define GUARD_SMEARING_CONTROL_H

#include <stdargs.h>

#include <buffers/gauge.h>
#include <buffers/adjoint.h>

#include <smearing/identity.h>
// #include <smearing/ape.h>
// #include <smearing/hyp.h>
#include <smearing/stout.h>
// #include <smearing/unitary.h>
// #include <smearing/hex.h>
// #include <smearing/gradient.h>


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
  smearing_type_e   type;

  /* Flags */
  int calculate_force_terms;
  int smearing_performed;
  
  /* Results -- main output for users */
  /* Both of these are shallow copies -- essentially pointers */
  gauge_field_t    result;
  adjoint_field_t  force_result;

  void* control;
} smearing_control;

smearing_control *construct_smearing_control(smearing_type type, ...);
void free_smearing_control(smearing_control *control);

#endif
