#include "utils.ih"

#include <global.h>
#include <start.h>
#include <sighandler.h>
#include <operator/tm_operators.h>
#include <linalg_eo.h>
#include <io/gauge.h>
#include <io/params.h>
#include <measure_gauge_action.h>
#include <ranlxd.h>
#include <read_input.h>
#include <expo.h>
#include <xchange/xchange.h>
#include <measure_rectangles.h>
#include <init/init_gauge_tmp.h>
#include <monomial/monomial.h>
#include <integrator.h>
#include <hamiltonian_field.h>
#include <update_tm.h>
#include <gettime.h>

#include <dirty_shameful_business.h>

void calculate_forces_numerically(su3adj *result, int const x, int const mu, int * mnllist, const int no)
{
  /* Set up the hamiltonian field */
  hamiltonian_field_t hf;
  hf.gaugefield = g_gauge_field;
  hf.momenta = 0;
  hf.derivative = 0;
  hf.update_gauge_copy = g_update_gauge_copy;
  hf.update_gauge_energy = g_update_gauge_energy;
  hf.update_rectangle_energy = g_update_rectangle_energy;
  hf.traj_counter = 0;
  
  /* Get some memory set aside for gauge fields and copy our current field */
  gauge_field_t original = get_gauge_field();
  memmove(original, g_gf, sizeof(su3_tuple) * (VOLUMEPLUSRAND + g_dbw2rand) + 1);
  
  gauge_field_t rotated = get_gauge_field();
   
  su3adj rotation;
  double *ar_rotation = (double*)&rotation;
  double *ar_result = (double*)result;
  double const epsilon = 5e-6;
   
  memmove(rotated, original, sizeof(su3_tuple) * (VOLUMEPLUSRAND + g_dbw2rand) + 1);
  
  su3 old_value;
  memmove(&old_value, &rotated[x][mu], sizeof(su3));
  for (int component = 0; component < 8; ++component)
  {
    /* Introduce a rotation along one of the components */
    memset(ar_rotation, 0, sizeof(su3adj));
    ar_rotation[component] = epsilon;
    
    su3 mat_rotation;
    exposu3(&mat_rotation, &rotation);
    
    _su3_times_su3(rotated[x][mu], mat_rotation, old_value);
    
    /* Calculate the action on the rotated field for the relevant monomials */
    double h_rotated_forward = 0.;
    
    for (int s_type = 0; s_type < no_smearings_monomial; ++s_type)
    {
      smear(smearing_control_monomial[s_type], rotated);
      ohnohack_remap_g_gauge_field(smearing_control_monomial[s_type]->result);
      
      for(int i = 0; i < no; ++i)
        if (monomial_list[ mnllist[i] ].smearing == s_type)
        {
          g_update_gauge_energy = 1;
          g_update_gauge_copy = 1;
          h_rotated_forward += monomial_list[ mnllist[i] ].accfunction(mnllist[i], &hf);
        }
    }
    ohnohack_remap_g_gauge_field(g_gf);
    
    memset(ar_rotation, 0, sizeof(su3adj));
    ar_rotation[component] = -epsilon;
    
    exposu3(&mat_rotation, &rotation);
    
    _su3_times_su3(rotated[x][mu], mat_rotation, old_value);
    
    /* Calculate the action on the rotated field for the relevant monomials */
    double h_rotated_backward = 0.;
    
    for (int s_type = 0; s_type < no_smearings_monomial; ++s_type)
    {
      smear(smearing_control_monomial[s_type], rotated);
      ohnohack_remap_g_gauge_field(smearing_control_monomial[s_type]->result);
      
      for(int i = 0; i < no; ++i)
        if (monomial_list[ mnllist[i] ].smearing == s_type)
        {
          g_update_gauge_energy = 1;
          g_update_gauge_copy = 1;
          h_rotated_backward += monomial_list[ mnllist[i] ].accfunction(mnllist[i], &hf);
        }
    }
    ohnohack_remap_g_gauge_field(g_gf);
    
    ar_result[component] = (h_rotated_forward - h_rotated_backward) / (2 * epsilon);
  }
  
  return_gauge_field(&rotated);
  return_gauge_field(&original);  
}
