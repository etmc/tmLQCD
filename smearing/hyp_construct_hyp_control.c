#include "hyp.ih"

hyp_control *construct_hyp_control(unsigned int iterations, double const coef_1, double const coef_2, double const coef_3)
{
  hyp_control *control = (hyp_control*)malloc(sizeof(hyp_control));
  control->coef[0] = coef_1;
  control->coef[1] = coef_2;
  control->coef[2] = coef_3;
  control->iterations = iterations;

  control->staples = (su3_two_index**)malloc(2 * sizeof(su3_two_index*));
  control->staples[0] = aalloc(VOLUMEPLUSRAND * sizeof(su3_two_index));
  control->staples[1] = aalloc(VOLUMEPLUSRAND * sizeof(su3_two_index));
  
  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  
  return control;
}
