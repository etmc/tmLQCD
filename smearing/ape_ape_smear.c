#include "ape.ih"

int ape_smear(gauge_field_t m_field_out, struct ape_parameters const *params, gauge_field_t m_field_in)
{
  /* We may alias the data, so we need something to store intermediate results somewhere else then m_field_out */
  gauge_field_t buffer = get_gauge_field();

  /* start of the the stout smearing **/
  for(int iter = 0; iter < params->iterations; ++iter)
  {
    generic_staples(buffer, m_field_in);
    APE_project_exclude_none(buffer, params->rho, buffer,  m_field_in);

    copy_gauge_field(m_field_out, buffer);
    exchange_gauge_field(m_field_out);

    m_field_in = m_field_out; /* Prepare for next iteration */
  }
  
  return_gauge_field(&buffer);

  return(0);
}
