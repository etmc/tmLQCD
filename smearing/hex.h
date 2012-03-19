#pragma once

#include <buffers/gauge.h>
#include <smearing/hyp.h>

/* Just to have a consistent look to the interface  */
typedef struct hyp_parameters hex_parameters;

/* All defined in terms of arrays of tuples -- needed to allow for g_gauge_field as input */
void stout_exclude_none(gauge_field_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
void stout_exclude_one (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
void stout_exclude_two (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);

int hex_smear(gauge_field_t m_field_out, hex_parameters const *params, gauge_field_t m_field_in);  /*  4 components in, 4 components out */
