#pragma once

#include <buffers/gauge.h>
#include <smearing/utils.h>

struct stout_parameters
{
  double rho;
  int    iterations;
};

void stout_smear(gauge_field_t m_field_out, struct stout_parameters const *params, gauge_field_t m_field_in);
void stout_smear_verbose(gauge_field_t m_field_out, struct stout_parameters const *params, gauge_field_t m_field_in);
void stout_links(gauge_field_t buff_out, double const coeff, gauge_field_t staples, gauge_field_t buff_in);
