#ifndef _TEST_CLOVER_SIX_INVERT_H
#define _TEST_CLOVER_SIX_INVERT_H

#include <cu/cu.h>

TEST(clover_six_invert);
TEST(clover_six_det);

TEST_SUITE(CLOVER){
  TEST_ADD(clover_six_invert),
    TEST_ADD(clover_six_det),
    TEST_SUITE_CLOSURE
};

#endif /* _TEST_CLOVER_SIX_INVERT_H */
