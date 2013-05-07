#ifndef GUARD_SMEARING_CONTROL_H
#define GUARD_SMEARING_CONTROL_H

#include <stdarg.h>

#include <buffers/gauge.h>
#include <buffers/adjoint.h>

#include <smearing/ape.h>
#include <smearing/ape_3d.h>
// #include <smearing/gradient.h>
#include <smearing/hex.h>
#include <smearing/hex_3d.h>
#include <smearing/hyp.h>
// #include <smearing/hyp_3d.h>
#include <smearing/identity.h>
// #include <smearing/jacobi.h>
#include <smearing/stout.h>
#include <smearing/stout_3d.h>
// #include <smearing/unitary.h>

typedef enum 
{
  APE,
  APE_3D,
  Gradient,
  HEX,
  HEX_3D,
  HYP,
  HYP_3D,
  Identity,
  Jacobi,
  Stout,
  Stout_3D,
  Unitary
} smearing_type;

extern char const * smearing_type_names[12];

typedef struct
{
    smearing_type type;
    int    id;
    int    iterations;
    double params[3]; 
    int    set[5];
} smearing_params_t;

typedef struct
{  
  smearing_type type;
  int id;
  
  /* Flags */
  int smearing_performed;
  
  /* Results -- main output for users */
  /* Both of these are shallow copies -- essentially just pointers */
  gauge_field_t    result;
  adjoint_field_t  force_result;

  void* type_control;
} smearing_control_t;

extern int no_smearings_declared;
extern smearing_params_t* smearing_declarations;

extern int no_smearings_monomial;
extern smearing_control_t **smearing_control_monomial;
  
extern int no_smearings_measurement;
extern smearing_control_t **smearing_control_measurement;

extern int no_smearings_operator;
extern smearing_control_t **smearing_control_operator;

smearing_control_t *construct_smearing_control(smearing_type type,  int calculate_force_terms, ...);
smearing_control_t *construct_smearing_control_from_params(smearing_params_t const *params,  int calculate_force_terms);
void free_smearing_control(smearing_control_t *control);

void smear(smearing_control_t *control, gauge_field_t in);
void smear_forces(smearing_control_t *control, adjoint_field_t in);

#endif
