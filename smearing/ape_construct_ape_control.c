#include "ape.ih"

ape_control *construct_ape_control(double rho, unsigned int iterations)
{
  ape_control *control = (ape_control*)malloc(sizeof(ape_control));
  control->rho = rho;
  control->iterations = iterations;
  control->smearing_performed = 0;

  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  return control;
}
