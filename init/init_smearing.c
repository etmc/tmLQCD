#include "init_smearing.h"
#include "fatal_error.h"

int no_smearings_monomial = 0;
smearing_control_t **smearing_control_monomial = NULL;
int no_smearings_measurement = 0;
smearing_control_t **smearing_control_measurement = NULL;
int no_smearings_operator = 0;
smearing_control_t **smearing_control_operator = NULL;

void init_smearing()
{
  /* The idea here is to convert the 'smearing declaration' array into a set of actual smearing control structs.
    * We can separate these out for the different applications. Monomials will need force calculations, for example,
    * whereas observables generally wouldn't. By setting up separate arrays containing just those smearing types that
    * are actually needed, we can later on just loop over these arrays without performing wasteful computations. The
    * freeing of the control structs is a secondary issue. But given how the information is now bound to monomials and
    * measurements, we can probably perform the clean up in those routines. After all, this would also guarantee that 
    * the smearing can be performed as long as the associated objects exist. */
  char error_string[256];
  
  for (int ctr = 0; ctr < no_monomials; ++ctr)
  {
    int id = monomial_list[ctr].smearing;
    int idx = 0;
    // Check if the smearing control already exists.
    while ((idx < no_smearings_monomial) && (id != smearing_control_monomial[idx]->id))
      ++idx;
    if (idx != no_smearings_monomial)
      monomial_list[ctr].smearing = idx;
    else
    {
      // If the smearing id does not yet exist, check if it has been declared at least
      idx = 0;
      while ((idx < no_smearings_declared) && (id != smearing_declarations[idx].id))
        ++idx;
      if (idx != no_smearings_declared)
      {
        smearing_control_t **old_array = smearing_control_monomial;
        smearing_control_monomial = (smearing_control_t**)calloc(no_smearings_monomial + 1, sizeof(smearing_control_t*));
        memmove(smearing_control_monomial, old_array, no_smearings_monomial * sizeof(smearing_control_t*));
        free(old_array);
        smearing_control_monomial[no_smearings_monomial] = construct_smearing_control_from_params(smearing_declarations + idx, 1 /* calculate_force_terms */);
        smearing_control_monomial[no_smearings_monomial]->id = id;
        monomial_list[ctr].smearing = no_smearings_monomial;
        ++no_smearings_monomial;
      }
      else
      {
        sprintf(error_string, "Smearing id %d given for monomial %d is not defined!", id, ctr);
        fatal_error(error_string, "cross_reference_smearing_ids");
      }
    }
  }
    
  for (int ctr = 0; ctr < no_measurements; ++ctr)
  {
    int id = measurement_list[ctr].smearing;
    int idx = 0;
    // Check if the smearing control already exists.
    while ((idx < no_smearings_measurement) && (id != smearing_control_measurement[idx]->id))
      ++idx;
    if (idx != no_smearings_measurement)
      measurement_list[ctr].smearing = idx;
    else
    {
      // If the smearing id does not yet exist, check if it has been declared at least
      idx = 0;
      while ((idx < no_smearings_declared) && (id != smearing_declarations[idx].id))
        ++idx;
      if (idx != no_smearings_declared)
      {
        smearing_control_t **old_array = smearing_control_measurement;
        smearing_control_measurement = (smearing_control_t**)calloc(no_smearings_measurement + 1, sizeof(smearing_control_t*));
        memmove(smearing_control_measurement, old_array, no_smearings_measurement * sizeof(smearing_control_t*));
        free(old_array);
        smearing_control_measurement[no_smearings_measurement] = construct_smearing_control_from_params(smearing_declarations + idx, 0 /* calculate_force_terms */);
        smearing_control_measurement[no_smearings_measurement]->id = id;
        measurement_list[ctr].smearing = no_smearings_measurement;
        ++no_smearings_measurement;
      }
      else
      {
        sprintf(error_string, "Smearing id %d given for measurement %d is not defined!", id, ctr);
        fatal_error(error_string, "cross_reference_smearing_ids");
      }
    }
  }
  
  for (int ctr = 0; ctr < no_operators; ++ctr)
  {
    int id = operator_list[ctr].smearing;
    int idx = 0;
    // Check if the smearing control already exists.
    while ((idx < no_smearings_operator) && (id != smearing_control_operator[idx]->id))
      ++idx;
    if (idx != no_smearings_operator)
      operator_list[ctr].smearing = idx;
    else
    {
      // If the smearing id does not yet exist, check if it has been declared at least
      idx = 0;
      while ((idx < no_smearings_declared) && (id != smearing_declarations[idx].id))
        ++idx;
      if (idx != no_smearings_declared)
      {
        smearing_control_t **old_array = smearing_control_operator;
        smearing_control_operator = (smearing_control_t**)calloc(no_smearings_operator + 1, sizeof(smearing_control_t*));
        memmove(smearing_control_operator, old_array, no_smearings_operator * sizeof(smearing_control_t*));
        free(old_array);
        smearing_control_operator[no_smearings_operator] = construct_smearing_control_from_params(smearing_declarations + idx, 0 /* calculate_force_terms */);
        smearing_control_operator[no_smearings_operator]->id = id;
        operator_list[ctr].smearing = no_smearings_operator;
        ++no_smearings_operator;
      }
      else
      {
        sprintf(error_string, "Smearing id %d given for operator %d is not defined!", id, ctr);
        fatal_error(error_string, "cross_reference_smearing_ids");
      }
    }
  }
  
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