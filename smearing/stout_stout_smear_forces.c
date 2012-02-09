#include "stout.ih"

void stout_smear_forces(struct stout_control *control, adjoint_field_t in)
{
  if (!control->calculate_force_terms)
    fatal_error("Stout control structure not setup for calculating force terms.", "stout_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Stout smearing not yet performed.", "stout_smear_forces");
  
  for (unsigned int ctr = 0; ctr < 3; ++ctr)
    control->scratch[ctr] = get_gauge_field();
  
  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
      _make_su3(control->scratch[0].field[x][mu], in.field[x][mu]);
  
  _Complex double trace;
  for (int iter = control->iterations - 1; iter > 0; --iter)
  {
    /* NOTE Write routines operating on whole gauge fields? This is clunky and tedious. */
    
    /* Calculate Gamma as defined by Peardon and Morningstar */
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        /* scratch[3] == Tr[Sigma' * B2 * U] * Q */
        _su3_times_su3(control->scratch[2].field[x][mu], control->scratch[0].field[x][mu], control->B2[iter].field[x][mu]);
        _trace_su3_times_su3(trace, control->scratch[2].field[x][mu], control->U[iter].field[x][mu]); /* NOTE Iteration of U before last! */
        _complex_times_su3(control->scratch[3].field[x][mu], trace, control->Q[iter].field[x][mu]);
        
        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q */
        _su3_times_su3(control->scratch[2].field[x][mu], control->scratch[0].field[x][mu], control->B1[iter].field[x][mu]);
        _trace_su3_times_su3(trace, control->scratch[2].field[x][mu], control->U[iter].field[x][mu]); /* NOTE Iteration of U before last! */
        _su3_add_equals_complex_identity(control->scratch[3].field[x][mu], trace)
        _su3_times_su3(control->scratch[1].field[x][mu], control->scratch[3].field[x][mu], control->Q[iter].field[x][mu]);

        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q  + f_1 * U * Sigma' */
        _su3_times_su3(control->scratch[2].field[x][mu], control->U[iter].field[x][mu], control->scratch[0].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], control->f1[iter][x][mu], control->scratch[2].field[x][mu]);
        _su3_plus_su3(control->scratch[1], control->scratch[3].field[x][mu], control->scratch[1]);
        
        /* scratch[1] == ( Tr[Sigma' * B1 * U] + Tr[Sigma' * B2 * U] * Q ) * Q  + f_1 * U * Sigma' + f2 * Q * U * Sigma' */
        _su3_times_su3(control->scratch[3].field[x][mu], control->Q[iter].field[x][mu], control->scratch[2].field[x][mu]);
        _complex_times_su3(control->scratch[3].field[x][mu], control->f2[iter][x][mu], control->scratch[3].field[x][mu]);
      }
  }
  
  for (unsigned int ctr = 0; ctr < 3; ++ctr)
    return_gauge_field(&control->scratch[ctr]);
}
