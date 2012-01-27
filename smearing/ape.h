#pragma once

#include <smearing/utils.h>

struct ape_parameters
{
  double rho;
  int    iterations;
};

int ape_smear(gauge_field_t m_field_out, struct ape_parameters const *params, gauge_field_t m_field_in);
void recombine_ape_smeared_tuples(gauge_field_t smeared, gauge_field_t original);
