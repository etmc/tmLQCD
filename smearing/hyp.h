#pragma once

#include <smearing/utils.h>

struct hyp_parameters
{
  double alpha[3];
  int    iterations;
};

/* All defined in terms of arrays of tuples -- needed to allow for g_gauge_field as input */

void hyp_staples_exclude_none(su3_tuple **buff_out, su3_tuple **buff_in); /* 12 components in, 12 components out */
void hyp_staples_exclude_one (su3_tuple **buff_out, su3_tuple **buff_in);  /* 12 components in, 12 components out */
void hyp_staples_exclude_two (su3_tuple **buff_out, su3_tuple  *buff_in);  /*  4 components in, 12 components out */

void APE_project_exclude_none(su3_tuple  *buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void APE_project_exclude_one (su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void APE_project_exclude_two (su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);

int hyp_smear(su3_tuple *m_field_out, struct hyp_parameters const *params, su3_tuple *m_field_in);  /*  4 components in, 4 components out */
