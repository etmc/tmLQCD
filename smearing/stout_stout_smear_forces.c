#include "stout.ih"

static void add_stout_terms_to_forces(gauge_field_t sigma, gauge_field_t const lambda, stout_notes_tuple const *trace, gauge_field_t const U)
{
  su3 staples;
  su3 aggr;
  
  for (unsigned int x = 0; x < VOLUME; ++x)
    for (unsigned int mu = 0; mu < 4; ++mu)
    {
      /* sigma' = (sigma' * expiQ) */
      aggr = sigma.field[x][mu];
      _su3_times_su3(sigma.field[x][mu], aggr, trace[x][mu].expiQ);
      
      /* sigma' = (sigma' * expiQ) + i * C^dag * lambda */
      generic_staples(&staples, x, mu, control->U[iter])
      _su3d_times_su3_acc(aggr, staples, lambda);
      _su3_imfac_acc(sigma.field[x][mu], 1, t1);

      /* sigma' = (sigma' * expiQ) + i * C^dag * lambda - i rho \sum [St.Deriv] */
      for (unsigned int nu = 0; nu < 4; ++nu)
        if (mu != nu)
          add_staple_derivatives(&aggr, mu, nu, lambda, U);
      _su3_imfac_acc(sigma.field[x][mu], -control->rho, aggr);
    }
}

void stout_smear_forces(stout_control *control, adjoint_field_t in)
{
  /* Check sanity of the call */
  if (!control->calculate_force_terms)
    fatal_error("Stout control structure not setup for calculating force terms.", "stout_smear_forces");

  if (!control->smearing_performed)
    fatal_error("Stout smearing not yet performed.", "stout_smear_forces");

  /* We'll need the forces in their tangent space representation, so let's first build this up. */
  gauge_field_t sigma = get_gauge_field();
  adjoint_to_fundamental(sigma, in);

  /* The staple derivatives contain non-local factors of lambda, so we'll need to store a whole field of those. */
  gauge_field_t lambda = get_gauge_field();

  
  /* The modifications are done backwards, all the time peeling off one layer of stouting... */
  for (unsigned sint iter = control->iterations - 1; iter > 0; --iter)
  {
    construct_lambda(lambda, sigma, control->trace[iter], control->U[iter]);
    add_stout_terms_to_forces(sigma, lambda, control->trace[iter], control->U[iter]);
  }

  /* The force terms are still in the tangent space representation, so project them back to the adjoint one */
  fundamental_to_adjoint(control->force_result, sigma);

  return_gauge_field(&lambda)
  return_gauge_field(&sigma);
}
