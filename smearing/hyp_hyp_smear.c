#include "hyp.ih"

#include "hyp_smear_first_stage.static"
#include "hyp_smear_second_stage.static"
#include "hyp_smear_third_stage.static"

void hyp_smear(hyp_control *control, gauge_field_t in)
{
  control->U[0] = in;
  
  for (unsigned int iter = 0; iter < control->iterations; ++iter)
  {
    smear_first_stage(control->staples[0], control->alpha[0], in);
    smear_second_stage(control->staples[1], control->alpha[1], control->staples[0], in);
    smear_third_stage(control->U[1], control->alpha[2], control->staples[1], in);
    in = control->U[1];
  }
  control->smearing_performed = 1;
}
