#include "stout.ih"

int stout_smear(su3_tuple *m_field_out, struct stout_parameters const *params, su3_tuple *m_field_in)
{
  static int initialized = 0;
  static su3_tuple *buffer;
  static su3 tmp;
  
  if (!initialized)
  {
    /* Allocate consecutive memory for both of the buffers upon first instantiation */
    buffer = (su3_tuple*)malloc(sizeof(su3_tuple) * VOLUMEPLUSRAND + 1);
#if (defined SSE || defined SSE2 || defined SSE3)
    buffer = (su3_tuple*)(((unsigned long int)(buffer) + ALIGN_BASE) & ~ALIGN_BASE);
#endif
    
    if (buffer == (su3_tuple*)NULL)
      return -1;
    initialized = 1;
  }

  /* start of the the stout smearing **/
  for(int iter = 0; iter < params->iterations; ++iter)
  {
    for (int x = 0; x < VOLUME; ++x)
      for (int mu = 0; mu < 4; ++mu)
      {
        generic_staples(&tmp, x, mu, m_field_in);
        _real_times_su3(tmp, params->rho, tmp);
        _su3_times_su3d(buffer[x][mu], tmp, m_field_in[x][mu]);
        project_antiherm(&buffer[x][mu]);
        exposu3_in_place(&buffer[x][mu]);
      }
    
    for(int x = 0; x < VOLUME; ++x)
      for(int mu = 0 ; mu < 4; ++mu)
      { 
        /* Input and output are allowed to be aliases -- use tmp */
        _su3_times_su3(tmp, buffer[x][mu], m_field_in[x][mu]);
        _su3_assign(m_field_out[x][mu], tmp);
      }

//    generic_exchange(m_field_out, sizeof(su3_tuple));
    m_field_in = m_field_out; /* Prepare for next iteration */
  }

  return(0);
}
