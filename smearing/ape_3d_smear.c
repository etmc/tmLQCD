#include "ape_3d.ih"

void ape_3d_smear(ape_3d_control *control, gauge_field_t in)
{
  /* We may alias the data, so we need something to store intermediate results somewhere else then m_field_out */
  control->U[0] = in;

  su3 staples;
  /* We need to take staples, so we need some working memory... */
  gauge_field_t buffer = get_gauge_field();

  double const rho_principal = 1.0 - control->rho;
  double const rho_staples = control->rho / 4.0;
  
  /* start of the the stout smearing **/
  for(unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    for (unsigned int x = 0; x < VOLUME; ++x)
    {
      _su3_assign(buffer[x][0], in[x][0]); // Left untouched, but still needed for future calculations!
      for (unsigned int mu = 1; mu < 4; ++mu)
      {
        generic_staples_3d(&staples, x, mu, in);
        _real_times_su3_plus_real_times_su3(buffer[x][mu], rho_principal, in[x][mu], rho_staples, staples);
        reunitarize(&buffer[x][mu]);
      }
    }

    /* Prepare for the next iteration -- the last result is now input! */
    swap_gauge_field(&control->U[1], &buffer);
    exchange_gauge_field(&control->U[1]);
    in = control->U[1];
  }
  control->result = control->U[1];
  return_gauge_field(&buffer);
}
