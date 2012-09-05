#ifndef _TEST_SU3_ALGEBRA_H
#define _TEST_SU3_ALGEBRA_H

#include <cu/cu.h>

TEST(su3_assign);
TEST(su3_expo_positivedet);
TEST(su3_multiply);
TEST(su3_inverse_multiply);
TEST(vector_add);
TEST(vector_sub);
TEST(vector_i_add);
TEST(vector_i_sub);
TEST(cmplx_times_vector);
TEST(cmplxcjg_times_vector);

TEST_SUITE(SU3_ALGEBRA){
  TEST_ADD(su3_assign),
    TEST_ADD(su3_expo_positivedet),
    TEST_ADD(su3_multiply),
    TEST_ADD(su3_inverse_multiply),
    TEST_ADD(vector_add),
    TEST_ADD(vector_sub),
    TEST_ADD(vector_i_add),
    TEST_ADD(vector_i_sub),
    TEST_ADD(cmplx_times_vector),
    TEST_ADD(cmplxcjg_times_vector),
    TEST_SUITE_CLOSURE
    };

#endif /* _TEST_SU3_ALGEBRA_H */


