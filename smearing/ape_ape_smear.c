#include "ape.ih"

int ape_smear(gauge_field_t m_field_out, struct ape_parameters const *params, gauge_field_t m_field_in)
{
  static su3_tuple tmp;
  double const rho_p = 1 - params->rho;
  double const rho_s = params->rho / 6.0;

  /* We may alias the data, so we need something to store intermediate results */
  gauge_field_t buffer = get_gauge_field();

  /* start of the the stout smearing **/
  for(int iter = 0; iter < params->iterations; ++iter)
  {
    generic_staples(buffer, m_field_in);
    recombine_ape_smeared_staples(buffer, m_field_in);

    copy_gauge_field(m_field_out, buffer);
    exchange_gauge_field(m_field_out);

    m_field_in = m_field_out; /* Prepare for next iteration */
  }
  
  return_gauge_field(&buffer);

  return(0);
}



