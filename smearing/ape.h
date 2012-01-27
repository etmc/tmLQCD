#pragma once

#include <smearing/utils.h>

struct ape_parameters
{
  double rho;
  int    iterations;
};

int ape_smear(gauge_field_t m_field_out, struct ape_parameters const *params, gauge_field_t m_field_in);
void APE_project_exclude_none(gauge_field_t buff_out, double const coeff, gauge_field_t staples, gauge_field_t buff_in);

