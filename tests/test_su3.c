
#include "test_su3_algebra.h"

TEST_SUITES {
  TEST_SUITE_ADD(SU3_ALGEBRA),
  TEST_SUITES_CLOSURE
};

int main(int argc,char *argv[]){
  CU_SET_OUT_PREFIX("regressions/");
  CU_RUN(argc,argv);
  return 0;
}


