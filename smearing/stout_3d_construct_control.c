#include "stout_3d.ih"

stout_3d_control *construct_stout_3d_control(int calculate_force_terms, unsigned int iterations, double rho)
{
  stout_3d_control *control = (stout_3d_control*)malloc(sizeof(stout_3d_control));
  control->rho = rho;
  control->iterations = iterations;
  
  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field_annotated("STOUT CONTROL U[1]");
  
  control->result = control->U[1]; /* A shallow copy, just putting the reference in place. */
  return control;
}
