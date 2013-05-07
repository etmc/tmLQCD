#include "ape.ih"

void ape_smear(ape_control *control, gauge_field_t in)
{
  /* We may alias the data, so we need something to store intermediate results somewhere else then m_field_out */
  control->U[0] = in;

  su3 staples;
  /* We need to take staples, so we need some working memory... */
  gauge_field_t buffer = get_gauge_field();

  double const coef_principal = 1.0 - control->coef;
  double const coef_staples = control->coef / 6.0;

  /* start of the the stout smearing **/
  for(unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        generic_staples(&staples, x, mu, in);
        _real_times_su3_plus_real_times_su3(buffer[x][mu], coef_principal, in[x][mu], coef_staples, staples);
        reunitarize(&buffer[x][mu]);
      }

    /* Prepare for the next iteration -- the last result is now input! */
    swap_gauge_field(&control->U[1], &buffer);
    exchange_gauge_field(&control->U[1]);
    in = control->U[1];
  }

  control->result = control->U[1];
  
  return_gauge_field(&buffer);
}
