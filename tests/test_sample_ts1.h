#ifndef _TEST_SAMPLE_TS1_H
#define _TEST_SAMPLE_TS1_H

#include <cu/cu.h>

/* test declarations, definitions in %.c */
TEST(test_true);
TEST(test_false);
TEST(test_fail);

/* define test suite (enum type thing) */
TEST_SUITE(TS1) {
  TEST_ADD(test_true),
  TEST_ADD(test_false),
  TEST_ADD(test_fail),
  TEST_SUITE_CLOSURE
};

#endif /* _TEST_SAMPLE_TS1_H */

