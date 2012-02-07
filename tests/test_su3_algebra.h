#ifndef _TEST_SU3_ALGEBRA_H
#define _TEST_SU3_ALGEBRA_H

#include <cu/cu.h>

TEST(su3_assign);
TEST(su3_expo_positivedet);

TEST_SUITE(SU3_ALGEBRA){
  TEST_ADD(su3_assign),
  TEST_ADD(su3_expo_positivedet),
  TEST_SUITE_CLOSURE
};

#endif /* _TEST_SU3_ALGEBRA_H */


