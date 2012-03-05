#include "hyp.ih"

void hyp_smear(hyp_control *control, gauge_field_t in)
{
  /* NOTE Can we get rid of one of these?  */
  gauge_field_t gamma[3];
  gauge_field_t v[3];

  for (unsigned int iter = 0; iter < params->iterations; ++iter)
  {
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        /* First level of contractions */
        hyp_staples_exclude_two(gamma, m_field_in);
        APE_project_exclude_two(v, params->alpha[2], gamma, m_field_in);
        exchange_gauge_field_array(v);

        /* Second level of contractions */
        hyp_staples_exclude_one(gamma, v);
        APE_project_exclude_one(v, params->alpha[1], gamma, m_field_in);
        exchange_gauge_field_array(v);

        /* Final level of contractions  */
        hyp_staples_exclude_none(gamma.field_array[0], v);
        APE_project_exclude_none(m_field_out, params->alpha[0], gamma.field_array[0], m_field_in);
        exchange_gauge_field(m_field_out);

        m_field_in = m_field_out; /* Prepare for next iteration */
      }
  }
  control->smearing_performed = 1;
  return_gauge_field_array(&v);
  return_gauge_field_array(&gamma);
  
  return 0;
}

