#include "stout.ih"

stout_control *construct_stout_control(double rho, unsigned int iterations, int calculate_force_terms)
{
  stout_control *control = (stout_control*)malloc(sizeof(stout_control));
  control->rho = rho;
  control->iterations = iterations;
  control->current = 0;
  control->calculate_force_terms = calculate_force_terms;

  control->scratch = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
  if (!calculate_force_terms)
  {
    control->U = (gauge_field_t*)malloc(2 * sizeof(gauge_field_t));
    control->U[1].field = get_gauge_field();
    result = control->U[1];
    return;
  }

  control->U  = (gauge_field_t*)malloc((iterations + 1) * sizeof(gauge_field_t));
  control->Q  = (gauge_field_t*)malloc(iterations * sizeof(gauge_field_t));
  control->B1 = (gauge_field_t*)malloc(iterations * sizeof(gauge_field_t));
  control->B2 = (gauge_field_t*)malloc(iterations * sizeof(gauge_field_t));
  
  control->f1 = (exp_par**)malloc(iterators * sizeof(exp_par*));
  control->f2 = (exp_par**)malloc(iterators * sizeof(exp_par*));
  
  for (unsigned int iter = 0; iter < iterations; ++iter)
  {
    control->U[iter + 1] = get_gauge_field();
    control->Q[iter]     = get_gauge_field();
    control->B1[iter]    = get_gauge_field();
    control->B2[iter]   = get_gauge_field();
    
    control->f1[iter] = (exp_par*)malloc(VOLUME * sizeof(exp_par));
    control->f2[iter] = (exp_par*)malloc(VOLUME * sizeof(exp_par));
  }
  
  result = control->U[iterations];
}
