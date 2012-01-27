#include "ape.ih"

int ape_smear(gauge_field_t m_field_out, struct ape_parameters const *params, gauge_field_t m_field_in)
{
  static su3_tuple tmp;
  double const rho_p = 1 - params->rho;
  double const rho_s = params->rho / 6.0;
  gauge_field_t buffer = get_gauge_field();

  /* start of the the stout smearing **/
  for(int iter = 0; iter < params->iterations; ++iter)
  {
    generic_staples(buffer, m_field_in);
    for (int x = 0; x < VOLUME; ++x)
      for (int mu = 0; mu < 4; ++mu)
      {
        _real_times_su3_plus_real_times_su3(buffer.field[x][mu], rho_p, m_field_in.field[x][mu], rho_s, tmp)
        reunitarize(&buffer.field[x][mu]);
      }
    
    for(int x = 0; x < VOLUME; ++x)
      for(int mu = 0 ; mu < 4; ++mu)
      {
        _su3_assign(m_field_out.field[x][mu], buffer.field[x][mu]);
      }

    generic_exchange(m_field_out.field, sizeof(su3_tuple));
    m_field_in = m_field_out; /* Prepare for next iteration */
  }
  
  return_gauge_field(&buffer);

  return(0);
}



