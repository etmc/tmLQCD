
#include "test_sample_ts1.h"
#include "test_sample_ts2.h"

TEST_SUITES {
  TEST_SUITE_ADD(TS1),
  TEST_SUITE_ADD(TS2),
  TEST_SUITES_CLOSURE
};

int main(int argc, char *argv[]) {
  CU_SET_OUT_PREFIX("regressions/");
  CU_RUN(argc,argv);

  return 0;
}


