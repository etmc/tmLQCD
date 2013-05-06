#include "hex_3d.ih"

#include "stout_fatten_links.static"
#include "hex_3d_smear_stage_1.static"
#include "hex_3d_smear_stage_2.static"

void hex_3d_smear(hex_3d_control *control, gauge_field_t in)
{
  control->U[0] = in; /* Shallow copy intended */

#pragma omp parallel
  for (int iter = 0; iter < control->iterations; ++iter)
  {
    smear_stage_1(control->V,    control->alpha[0],             in);
    smear_stage_2(control->U[1], control->alpha[1], control->V, in);

#pragma omp single
    in = control->U[1]; /* Shallow copy intended */
  }
  
  control->result = control->U[1];
}
