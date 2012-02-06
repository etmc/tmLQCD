#ifndef _TEST_SAMPLE_TS2_H
#define _TEST_SAMPLE_TS2_H

#include <cu/cu.h>

/* test declarations, definitions in %.c */
TEST(test_true2);

/* define test suite (enum type thing) */
TEST_SUITE(TS2) {
  TEST_ADD(test_true2),
  TEST_SUITE_CLOSURE
};

#endif /* _TEST_SAMPLE_TS2_H */

