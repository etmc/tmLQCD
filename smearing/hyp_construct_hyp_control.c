#include "hyp.ih"

hyp_control *construct_hyp_control(double const alpha_0, double const alpha_1, double const alpha_2,  unsigned int iterations)
{
  hyp_control *control = (ape_control*)malloc(sizeof(hyp_control));
  control->alpha[0] = alpha_0;
  control->alpha[1] = alpha_1;
  control->alpha[2] = alpha_2;
  control->iterations = iterations;
  control->smearing_performed = 0;

  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  return control;
}
