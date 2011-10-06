#include "hex.ih"

int hex_smear(su3_tuple *m_field_out, hex_parameters const *params, su3_tuple *m_field_in)
{
  static int initialized = 0;
  static su3_tuple *gamma_buffer[3];
  static su3_tuple *v_buffer[3];

  if (!initialized)
  {
    /* Allocate consecutive memory for both of the buffers upon first instantiation */
    /* Three times 4 buffers needed for compatibility purposes (similar signature to gauge_field...) */
    for (int idx = 0; idx < 3; ++idx)
    {
      gamma_buffer[idx] = (su3_tuple*)malloc(sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
      v_buffer[idx] = (su3_tuple*)malloc(sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
      if ((gamma_buffer[idx] == (su3_tuple*)NULL) || (v_buffer[idx] == (su3_tuple*)NULL))
        return -1;
#if (defined SSE || defined SSE2 || defined SSE3)
      gamma_buffer[idx] = (su3_tuple*)(((unsigned long int)(gamma_buffer[idx]) + ALIGN_BASE) & ~ALIGN_BASE);
      v_buffer[idx] = (su3_tuple*)(((unsigned long int)(v_buffer[idx]) + ALIGN_BASE) & ~ALIGN_BASE);
#endif
    }
    initialized = 1;
  }

  for (int iter = 0; iter < params->iterations; ++iter)
  {
    /* First level of contractions */
    hyp_staples_exclude_two(gamma_buffer, m_field_in);
    stout_exclude_two(v_buffer, params->alpha[2], gamma_buffer, m_field_in);
    for (int idx = 0; idx < 3; ++idx)
      generic_exchange(v_buffer[idx], sizeof(su3_tuple));

    /* Second level of contractions */
    hyp_staples_exclude_one(gamma_buffer, v_buffer);
    stout_exclude_one(v_buffer, params->alpha[1], gamma_buffer, m_field_in);
    for (int idx = 0; idx < 3; ++idx)
      generic_exchange(v_buffer[idx], sizeof(su3_tuple));

    /* Final level of contractions  */
    hyp_staples_exclude_none(gamma_buffer, v_buffer);
    stout_exclude_none(m_field_out, params->alpha[0], gamma_buffer, m_field_in);
    generic_exchange(m_field_out, sizeof(su3_tuple));
    
    m_field_in = m_field_out; /* Prepare for next iteration */
  }

  return 0;
}

