#pragma once

#include <buffers/gauge.h>
#include <smearing/hyp.h>

/* Just to have a consistent look to the interface  */
typedef struct hyp_parameters hex_parameters;

/* NOTE(gkanwar): Temporarily disable the new gauge_field_t interfaces to enable
 * compilation with old su3_tuple array signatures. */
#if 0
/* All defined in terms of arrays of tuples -- needed to allow for g_gauge_field as input */
void stout_exclude_none(gauge_field_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
void stout_exclude_one (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);
void stout_exclude_two (gauge_field_array_t buff_out, double const coeff, gauge_field_array_t staples, gauge_field_t buff_in);

int hex_smear(gauge_field_t m_field_out, hex_parameters const *params, gauge_field_t m_field_in);  /*  4 components in, 4 components out */
#endif


void stout_exclude_one(su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void stout_exclude_two(su3_tuple **buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);
void stout_exclude_none(su3_tuple *buff_out, double const coeff, su3_tuple **staples, su3_tuple *buff_in);

int hex_smear(su3_tuple *m_field_out, hex_parameters const *params, su3_tuple *m_field_in);
