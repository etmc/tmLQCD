#pragma once

#include <smearing/hyp.h>

/* Just to have a consistent look to the interface  */
typedef struct hyp_parameters hex_parameters;

/* All defined in terms of arrays of tuples -- needed to allow for g_gauge_field as input */
void stout_exclude_none(su3_tuple  *buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void stout_exclude_one (su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void stout_exclude_two (su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);

int hex_smear(su3_tuple *m_field_out, hex_parameters const *params, su3_tuple *m_field_in);  /*  4 components in, 4 components out */
