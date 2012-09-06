#define MAIN_PROGRAM

#include "test_linalg_spinor.h"

TEST_SUITES {
  TEST_SUITE_ADD(LINALG),
  TEST_SUITES_CLOSURE
};

int main(int argc,char *argv[]){
  CU_SET_OUT_PREFIX("regressions/");
  CU_RUN(argc,argv);
  return 0;
}

#undef MAIN_PROGRAM
