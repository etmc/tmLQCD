#include "hex_3d.ih"

hex_3d_control *construct_hex_3d_control(int calculate_force_terms, unsigned int iterations, double const alpha_1, double const alpha_2, double const alpha_3)
{
  hex_3d_control *control = (hex_3d_control*)malloc(sizeof(hex_3d_control));
  control->alpha[0] = alpha_1;
  control->alpha[1] = alpha_2;
  control->iterations = iterations;
  
  /* We can keep this quite simple if we don't need any forces anyway. */
  control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  control->U[1] = get_gauge_field();
  
  control->V = (su3_two_index  *)aalloc(sizeof(su3_two_index) * (VOLUMEPLUSRAND + g_dbw2rand));
  
  return control;
}
