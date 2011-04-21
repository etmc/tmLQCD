#pragma once

#include <smearing/utils.h>

struct ape_parameters
{
  double rho;
  int    iterations;
};

int ape_smear(su3_tuple *m_field_out, struct ape_parameters const *params, su3_tuple *m_field_in);