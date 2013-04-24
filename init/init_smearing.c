#include "init_smearing.h"
#include "fatal_error.h"

typedef enum
{
  FLAG_MONOMIAL,
  FLAG_MEASUREMENT,
  FLAG_OPERATOR
} type_flag_t;

int no_smearings_monomial = 0;
smearing_control_t **smearing_control_monomial = NULL;
int no_smearings_measurement = 0;
smearing_control_t **smearing_control_measurement = NULL;
int no_smearings_operator = 0;
smearing_control_t **smearing_control_operator = NULL;

static void build_targeted_array(type_flag_t type)
{
  char error_string[256];
  
  int *no_smearings;
  smearing_control_t ***smearing_control;
  int calculate_force_terms = 0;
  int *id = 0;
  
  switch (type)
  {
    case FLAG_MONOMIAL:
      no_smearings = &no_smearings_monomial;
      smearing_control = &smearing_control_monomial;
      calculate_force_terms = 1;
      break;
    case FLAG_MEASUREMENT:
      no_smearings = &no_smearings_measurement;
      smearing_control = &smearing_control_measurement;
      break;
    case FLAG_OPERATOR:
      no_smearings = &no_smearings_operator;
      smearing_control = &smearing_control_operator;
  };
  
  for (int ctr = 0; ctr < no_monomials; ++ctr)
  {
    switch (type)
    {
      case FLAG_MONOMIAL:    id = &(monomial_list[ctr].smearing);    break;
      case FLAG_MEASUREMENT: id = &(measurement_list[ctr].smearing); break;
      case FLAG_OPERATOR:    id = &(operator_list[ctr].smearing);
    }

    // Check if the smearing control already exists.
    int idx = 0;
    while ((idx < *no_smearings) && (*id != (*smearing_control)[idx]->id))
      ++idx;
    if (idx != *no_smearings)
      *id = idx;
    else
    {
      // If the smearing id does not yet exist, check if it has been declared at least
      idx = 0;
      while ((idx < no_smearings_declared) && (*id != smearing_declarations[idx].id))
        ++idx;
      if (idx != no_smearings_declared)
      {
        smearing_control_t **old_array = *smearing_control;
        *smearing_control = (smearing_control_t**)calloc(*no_smearings + 1, sizeof(smearing_control_t*));
        memmove((*smearing_control), old_array, *no_smearings * sizeof(smearing_control_t*));
        free(old_array);
        (*smearing_control)[*no_smearings] = construct_smearing_control_from_params(smearing_declarations + idx, calculate_force_terms);
        (*smearing_control)[*no_smearings]->id = *id;
        *id = *no_smearings;
        ++(*no_smearings);
      }
      else
      {
        sprintf(error_string, "Smearing id %d is not defined!", *id);
        fatal_error(error_string, "cross_reference_smearing_ids");
      }
    }
  }
}

void init_smearing()
{
  /* The idea here is to convert the 'smearing declaration' array into a set of actual smearing control structs.
    * We can separate these out for the different applications. Monomials will need force calculations, for example,
    * whereas observables generally wouldn't. By setting up separate arrays containing just those smearing types that
    * are actually needed, we can later on just loop over these arrays without performing wasteful computations. The
    * freeing of the control structs is a secondary issue. But given how the information is now bound to monomials and
    * measurements, we can probably perform the clean up in those routines. After all, this would also guarantee that 
    * the smearing can be performed as long as the associated objects exist. */
  build_targeted_array(FLAG_MONOMIAL);
  build_targeted_array(FLAG_MEASUREMENT);
  build_targeted_array(FLAG_OPERATOR);
  
  free(smearing_declarations);
  no_smearings_declared = 0;
}

void finalize_smearing()
{
  for (int ctr = 0; ctr < no_smearings_monomial; ++ctr)
    free_smearing_control(smearing_control_monomial[ctr]);
  free(smearing_control_monomial);
  for (int ctr = 0; ctr < no_smearings_measurement; ++ctr)
    free_smearing_control(smearing_control_measurement[ctr]);
  free(smearing_control_measurement);
  for (int ctr = 0; ctr < no_smearings_operator; ++ctr)
    free_smearing_control(smearing_control_operator[ctr]);
  free(smearing_control_operator);  
}