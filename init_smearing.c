#include <global.h>
#include <init_smearing.h>

smearing_control_t **smearing_control = NULL;
int no_smearing_types = 2; /* FIXME Hardcoded, but a future argument. */

void init_smearing(int no_smearing_types)
{
  /* FIXME We'll use a hardcoded setup for now, but anticipate multiple smearing constructs. */
  smearing_control = malloc(no_smearing_types * sizeof(smearing_control_t*)); 

  int ohnohack_stout_calculate_force = 1;
  double ohnohack_stout_rho = 0.9 / 6.0; /* Think about changing this to alpha -- we need to choose a convention. */
  unsigned int ohnohack_stout_no_iter = 6;

  smearing_control[0] = construct_smearing_control(Identity, ohnohack_stout_calculate_force);
  smearing_control[1] = construct_smearing_control(Stout, ohnohack_stout_calculate_force, ohnohack_stout_rho, ohnohack_stout_no_iter);
}

void finalize_smearing()
{
  /* FIXME This should be promoted to a flexible loop later on. */
  for (int ctr = 0; ctr < 2; ++ctr)
    free_smearing_control(smearing_control[ctr]);
  free(smearing_control);
}
