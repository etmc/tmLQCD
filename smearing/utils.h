#pragma once

#include <su3.h>
#include <buffers/gauge.h>

/* We need a number of indices to do the bookkeeping.
   This will always amount to at most 12 fields, but
   we define some aliases so that we don't have to do
   the mental mapping all the time. */

enum I0
{
  I0_0 = 0, I0_1 = 1, I0_2 = 2, I0_3 = 3
};

enum I1
{
  I1_0_1 =  0, I1_0_2 =  1, I1_0_3 =  2, I1_1_0 =  3,
  I1_1_2 =  4, I1_1_3 =  5, I1_2_0 =  6, I1_2_1 =  7,
  I1_2_3 =  8, I1_3_0 =  9, I1_3_1 = 10, I1_3_2 = 11
};

enum I2
{
  I2_0_12 =  0, I2_0_21 =  0, I2_0_13 =  1, I2_0_31 =  1,
  I2_0_23 =  2, I2_0_32 =  2, I2_1_02 =  3, I2_1_20 =  3,
  I2_1_03 =  4, I2_1_30 =  4, I2_1_23 =  5, I2_1_32 =  5,
  I2_2_01 =  6, I2_2_10 =  6, I2_2_03 =  7, I2_2_30 =  7,
  I2_2_13 =  8, I2_2_31 =  8, I2_3_01 =  9, I2_3_10 =  9,
  I2_3_02 = 10, I2_3_20 = 10, I2_3_12 = 11, I2_3_21 = 11
};

void generic_staples(gauge_field_t buff_out, int x, int mu, gauge_field_t buff_in);
void generic_exchange(void *field_in, int bytes_per_site);
void project_antiherm(su3 *omega);
void project_herm(su3 *omega);
void reunitarize(su3 *omega);

void print_su3(su3 *in);
void print_config_to_screen(su3 **in);
