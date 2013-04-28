#include "ape_3d.ih"

ape_3d_control *construct_ape_3d_control(unsigned int iterations, double rho)
{
  ape_3d_control *control = (ape_3d_control*)malloc(sizeof(ape_3d_control));
  control->rho = rho;
  control->iterations = iterations;
  control->smearing_performed = 0;

  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  return control;
}
