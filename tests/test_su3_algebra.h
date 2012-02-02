#include <cu/cu.h>

TEST(su3_assign);
TEST(su3_expo_positivedet);

TEST_SUITE(SU3_ALGEBRA){
  TEST_ADD(su3_assign),
  TEST_ADD(su3_expo_positivedet),
  TEST_SUITE_CLOSURE
};

