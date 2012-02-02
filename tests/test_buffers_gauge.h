#include <cu/cu.h>

TEST(buffers_gauge_allocate_finalize);
TEST(buffers_gauge_get_return);

TEST_SUITE(BUFFERS_GAUGE){
  TEST_ADD(buffers_gauge_allocate_finalize),
  TEST_ADD(buffers_gauge_get_return),
  TEST_SUITE_CLOSURE
};

