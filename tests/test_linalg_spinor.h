#ifndef _TEST_LINALG_SPINOR_H
#define _TEST_LINALG_SPINOR_H

#include <cu/cu.h>

TEST(scalar_prod_real);
TEST(snorm);
TEST(sdiff);
TEST(aaddm_r);
TEST(amadd_r);

TEST_SUITE(LINALG){
  TEST_ADD(scalar_prod_real),
    TEST_ADD(snorm),
    TEST_ADD(sdiff),
    TEST_ADD(aaddm_r),
    TEST_ADD(amadd_r),
    TEST_SUITE_CLOSURE
    };

#endif
