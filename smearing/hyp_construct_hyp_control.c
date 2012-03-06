#include "hyp.ih"

hyp_control *construct_hyp_control(double const alpha_1, double const alpha_2, double const alpha_3,  unsigned int iterations)
{
  hyp_control *control = (hyp_control*)malloc(sizeof(hyp_control));
  control->alpha[0] = alpha_1;
  control->alpha[1] = alpha_2;
  control->alpha[2] = alpha_3;
  control->iterations = iterations;
  control->smearing_performed = 0;

  control->staples = (su3_outer**)malloc(2 * sizeof(su3_outer*));
  control->staples[0] = aalloc(VOLUMEPLUSRAND * sizeof(su3_outer));
  control->staples[1] = aalloc(VOLUMEPLUSRAND * sizeof(su3_outer));
  
  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  
  return control;
}
