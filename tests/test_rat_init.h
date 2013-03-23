#ifndef _TEST_RAT_INIT_H
#define _TEST_RAT_INIT_H

#include <cu/cu.h>

TEST(rat_init);

TEST_SUITE(RAT){
  TEST_ADD(rat_init),
    TEST_SUITE_CLOSURE
};

#endif
