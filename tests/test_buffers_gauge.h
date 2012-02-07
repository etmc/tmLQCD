#ifndef _TEST_BUFFERS_GAUGE_H
#define _TEST_BUFFERS_GAUGE_H

#include <cu/cu.h>

TEST(buffers_gauge_allocate_finalize);
TEST(buffers_gauge_get_return);

TEST_SUITE(BUFFERS_GAUGE){
  TEST_ADD(buffers_gauge_allocate_finalize),
  TEST_ADD(buffers_gauge_get_return),
  TEST_SUITE_CLOSURE
};

#endif /* _TEST_BUFFERS_GAUGE_H */


