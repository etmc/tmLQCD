#pragma once

#include <smearing/utils.h>

struct stout_parameters
{
  double rho;
  int    iterations;
};

int stout_smear(gauge_field_t m_field_out, struct stout_parameters const *params, gauge_field_t m_field_in);