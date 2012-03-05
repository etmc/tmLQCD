#include "hyp.ih"

void hyp_smear(hyp_control *control, gauge_field_t in)
{
  control->U[0] = in;
  
  gauge_field_t buffer_1[3];
  gauge_field_t buffer_2[3];
  for (unsigned int ctr = 0; ctr < 3; ++ctr)
    gamma[ctr] = get_gauge_field();

  for (unsigned int iter = 0; iter < params->iterations; ++iter)
  {
    for (unsigned int x = 0; x < VOLUME; ++x)
      for (unsigned int mu = 0; mu < 4; ++mu)
      {
        hyp_smear_first_stage(buffer_1, control->alpha[0], in);
        hyp_smear_second_stage(buffer_2, control->alpha[1], buffer_1);
        hyp_smear_third_stage(control->U[1], control->alpha[2], buffer_2);
      }
    in = control->U[1];
  }
  control->smearing_performed = 1;

  for (unsigned int ctr = 0; ctr < 3; ++ctr)
  {
    return_gauge_field_array(&gamma[ctr]);
    return_gauge_field_array(&v[ctr]);
  }
  
  return 0;
}

