#include <cu/cu.h>

TEST(test_true) {
  assertTrue(1);
}

TEST(test_false) {
  assertFalse(0);
}

TEST(test_fail) {
  assertFalse(1);
}

