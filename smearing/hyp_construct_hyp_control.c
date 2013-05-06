#include "hyp.ih"

hyp_control *construct_hyp_control(unsigned int iterations, double const alpha_1, double const alpha_2, double const alpha_3)
{
  hyp_control *control = (hyp_control*)malloc(sizeof(hyp_control));
  control->alpha[0] = alpha_1;
  control->alpha[1] = alpha_2;
  control->alpha[2] = alpha_3;
  control->iterations = iterations;

  control->staples = (su3_two_index**)malloc(2 * sizeof(su3_two_index*));
  control->staples[0] = aalloc(VOLUMEPLUSRAND * sizeof(su3_two_index));
  control->staples[1] = aalloc(VOLUMEPLUSRAND * sizeof(su3_two_index));
  
  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  
  return control;
}
