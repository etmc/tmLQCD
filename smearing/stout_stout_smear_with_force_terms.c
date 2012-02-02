#include "stout.ih"

stout_force_intermediate *stout_smear_with_force_terms(gauge_field_t m_field_out, struct stout_parameters const *params, gauge_field_t m_field_in)
{
  gauge_field_t buffer = get_gauge_field();
  
  /* Initialize a buffer for the intermediate results.*/
  /* NOTE Look into avoiding this malloc/free cycle. The problem is: how to avoid
          memory leaks from unreturned fields. Maybe static parameters somehow? 
          How about subsequent calls with different numbers of iterations? */
  /* NOTE This should probably be separated in a function that creates and destroys the structure.
          We can then take it as a parameter to this function and assume its proper initialization. 
          This should give the optimal balance between user control and security. */
  stout_force_intermediate *m_force = (stout_force_intermediate *)malloc(sizeof(stout_force_intermediate));
  m_force->U  = (gauge_field_t*)malloc(params->iterations * sizeof(gauge_field_t));
  m_force->Q  = (gauge_field_t*)malloc(params->iterations * sizeof(gauge_field_t));
  m_force->f0 = (complex_field_array_t*)malloc(params->iterations * sizeof(complex_field_array_t));
  m_force->f1 = (complex_field_array_t*)malloc(params->iterations * sizeof(complex_field_array_t));
  m_force->f2 = (complex_field_array_t*)malloc(params->iterations * sizeof(complex_field_array_t));
  m_force->u  = (real_field_array_t*)malloc(params->iterations * sizeof(real_field_array_t));
  m_force->v  = (real_field_array_t*)malloc(params->iterations * sizeof(real_field_array_t));
  
  for (unsigned int iter = 0; iter < params->iterations; ++iter)
  {
    m_force->U[iter]  = get_gauge_field(); /* Will store the smeared field */
    m_force->Q[iter]  = get_gauge_field(); /* For the 'weighted plaquette' terms */
    m_force->f0[iter] = get_complex_field_array(4);
    m_force->f1[iter] = get_complex_field_array(4);
    m_force->f2[iter] = get_complex_field_array(4);
    m_force->u[iter]  = get_real_fsield_array(4);
    m_force->v[iter]  = get_real_field_array(4);
  }
  
  /* Start of the the stout smearing */
  for(int iter = 0; iter < params->iterations; ++iter)
  {
    generic_staples(m_force->Q[iter], m_field_in);
    stout_links_with_force_terms(buffer, m_force->Q[iter], m_force->f0[iter], m_force->f1[iter], m_force->f2[iter],
                                         m_force->u[iter], m_force->v[iter], params->rho, buffer, buff_in)
    copy_gauge_field(m_field_out, buffer);
    exchange_gauge_field(m_field_out);

    m_field_in = m_field_out; /* Prepare for next iteration */
  }
  
  return_gauge_field(&buffer);
  
  return m_force;
}
